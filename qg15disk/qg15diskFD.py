"""
Solves the nondimensional 1.5-layer QG equations on a gamma plane

Q = lap(psi) - 1/Bu*psi - 1/Rog*r^2
D/Dt(Q)=0

Burger number Bu = (Ld/L)^2 
Curvature rossby number Rog = 2U/(gamma*L^3) = 2U/(fL)*(ap^2/L^2) = inertial/curvature

The size of the polar cap is also important to segregate cyclones and anticyclones

To run and plot using e.g. 8 processes:
    
    RUN="t1"
    DIR="FDoutput"
    mpiexec -n 8 python3 qg15diskFD.py $RUN $DIR
    
    python3 plot_energeticsFD.py $RUN $DIR

    rm -rf frames
    mpiexec -n 8 python3 plot_diskFD.py $RUN $DIR snapshots/*.h5
    ffmpeg -y -framerate 24 -i frames/write_%06d.png -c:v libx264 -pix_fmt yuv420p $DIR/${RUN}_output.mp4

"""

import csv
import pathlib
import sys
import numpy as np
import dedalus.public as d3
from dedalus.core.domain import *
import os, re, glob
from glob import glob
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # safe on headless clusters
import matplotlib.pyplot as plt
import logging
logger = logging.getLogger(__name__)
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.rank

# declare runname
run_name = sys.argv[1] if len(sys.argv) > 1 else "default"
outdir   = sys.argv[2] if len(sys.argv) > 2 else f"output_{run_name}"
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

def safe_scalar(expr):
    val = expr.evaluate().allgather_data()
    return val.item() if val.size > 0 else 0.0

def extract_stepnum(fname):
    m = re.search(r'_s(\d+)', fname)
    return int(m.group(1)) if m else -1

# physical parameters
Bu = 1
Rog = 10
nu = 1e-9
mu = 1e-2

# Simulation parameters
R = 5 # size of domain relative to physical length scale L
Nphi = 256
Nr = int(np.round(Nphi/np.pi/4)*4)
dealias = 3/2
stop_sim_time = 500
write_speed = stop_sim_time/500
timestepper = d3.SBDF1
max_timestep = 1
dtype = np.float64
    # want Nphi at least 2*pi*R*N/sqrt(Bu) but a power of 2 for performance
    # number of cells per deformation radius N approx 10 if possible
    # Nr approx Nphi/pi but rounded to a multiple of four

# Bases
coords = d3.PolarCoordinates('phi', 'r')
dist = d3.Distributor(coords, dtype=dtype) # change mesh based on number of processes
disk = d3.DiskBasis(coords, shape=(Nphi, Nr), radius=R, dealias=dealias, dtype=dtype)
edge = disk.edge
domain = Domain(dist, bases=[disk])

# Fields
t = dist.Field()
deltaT = dist.Field()
iterNum = dist.Field()
eta = dist.Field(name='eta', bases=disk)
noise = dist.Field(name='noise', bases=disk)

psi = dist.Field(name='psi', bases=disk)
q0 = dist.Field(name='q0',bases=disk)
qf = dist.Field(name='qf', bases=disk)
qOU = dist.Field(name='qOU', bases=disk)

q = dist.Field(name='q', bases=disk)
a1 = dist.Field(name='a1', bases=disk)   # a1 = lap(q)
a2 = dist.Field(name='a2', bases=disk)   # a2 = lap(a1)
a3 = dist.Field(name='a3', bases=disk)   # a3 = lap(a2)

tau_psi = dist.Field(name='tau_psi', bases=edge)
tau_q  = dist.Field(name='tau_q',  bases=edge)
tau_a1 = dist.Field(name='tau_a1', bases=edge)
tau_a2 = dist.Field(name='tau_a2', bases=edge)
tau_a3 = dist.Field(name='tau_a3', bases=edge)

# Substitutions
phi, r = dist.local_grids(disk)
lift = lambda A: d3.Lift(A, disk, -1)
dr = lambda A: d3.radial(d3.grad(A)(r=R))
u = d3.Skew(d3.Gradient(psi))
eta['g'] = -1/Rog*r**2

# hyperdiffusion (4th order) and drag (forms following Gallet and Ferrari 2021) 
HD = nu*d3.lap(a3)
drag = - mu*d3.div(np.sqrt(d3.grad(psi)@d3.grad(psi))*d3.grad(psi))
    # make sure drag arrest wavenumber is significantly smaller (at least 2 orders) than kd=2*pi*L/Ld

# build helmholtz filter for noise
ell = 4*np.pi*R/Nphi # forcing scale in physical space (at 2x Nyquist)
problem = d3.LBVP([qf, a1, a2, a3, tau_q, tau_a1, tau_a2, tau_a3],namespace=locals())
problem.add_equation("a1 - lap(qf)  + lift(tau_a1) = 0")
problem.add_equation("a2 - lap(a1)  + lift(tau_a2) = 0")
problem.add_equation("a3 - lap(a2)  + lift(tau_a3) = 0")
problem.add_equation("qf - ell**8*a3 + lift(tau_q) = q0")
problem.add_equation("qf(r=R) = 0")
problem.add_equation("a1(r=R) = 0")
problem.add_equation("a2(r=R) = 0")
problem.add_equation("a3(r=R) = 0")
forcing_filter = problem.build_solver()


# FORCING FUNCTION
target_rms = 1
q0.change_scales(dealias)
qf.change_scales(dealias)
theta = 1  # OU correlation time
def stochastic_forcing(*args):
    deltaT = args[1].allgather_data().item()
    iterNum = args[2].allgather_data().item()
    rho = np.exp(-deltaT/theta)
    
    # fill forcing field
    seedVal=int(42+iterNum*7919) # deterministic
    q0.fill_random('g', seed=seedVal, distribution='normal', scale=1)
    forcing_filter.solve()

    # rms rescale
    rms = np.sqrt(safe_scalar(d3.integ(qf*qf)))
    qf['g'] *= np.sqrt(target_rms**2 / rms**2)

    # OU update
    qOU.change_scales(dealias)
    qOU['g'] = rho*qOU['g'] + np.sqrt(max(1-rho**2,0))*qf['g']
    
    return qOU['g']

def S(*args, domain=domain, F=stochastic_forcing):
    return d3.GeneralFunction(
        dist=dist,
        domain=domain,
        tensorsig=(),
        dtype=np.float64,
        layout="g",
        func=F,
        args=args,
    )

forcing = S(t,deltaT,iterNum)

# Main problem
problem = d3.IVP([q, psi, a1, a2, a3, tau_q, tau_a1, tau_a2, tau_a3, tau_psi],time=t, namespace=locals())
problem.add_equation("dt(q) + HD + lift(tau_q) = - u@grad(q + eta) + forcing + drag")
problem.add_equation("a1 - lap(q)  + lift(tau_a1) = 0")
problem.add_equation("a2 - lap(a1) + lift(tau_a2) = 0")
problem.add_equation("a3 - lap(a2) + lift(tau_a3) = 0")
problem.add_equation("q - lap(psi) + 1/Bu*psi + lift(tau_psi) = 0")

problem.add_equation("psi(r=R) = 0") # material surface
problem.add_equation("q(r=R)  = 0") # free slip
problem.add_equation("a1(r=R)  = 0")
problem.add_equation("a2(r=R)  = 0")
problem.add_equation("a3(r=R)  = 0")

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time
file_handler_mode = 'overwrite'
initial_timestep = max_timestep/10

# Analysis
Qtot=q+eta # total PV
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=write_speed, max_writes=10, mode=file_handler_mode)
snapshots.add_task(q, name='q')
snapshots.add_task(d3.lap(psi), name='zeta')
snapshots.add_task(Qtot, name='Qtot')
snapshots.add_task(psi, name='psi')

# CFL
CFL = d3.CFL(solver, initial_dt=initial_timestep, cadence=10, safety=0.5, threshold=0.1,
             max_change=1.5, min_change=0.5, max_dt=max_timestep)
CFL.add_velocity(u)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=10)
flow.add_property((psi*psi), name='psi2')
flow.add_property((u@u), name='u2')
flow.add_property((q*q), name='q2')

# write parameters to CSV
params = {
    "Bu": Bu,
    "Rog": Rog,
    "nu": nu,
    "mu": mu,
    "R": R,
    "Nr": Nr,
    "Nphi": Nphi,
    "forcing target rms": target_rms,
    "forcing OU correlation time": theta,
    "forcing filter scale": ell,
}

if rank == 0:
    with open(f"{outdir}/{run_name}_parameters.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["parameter", "value"])
        for key, val in params.items():
            writer.writerow([key, val])


# prepare diagnostics CSV
diagnostics_csv = os.path.join(outdir, f"{outdir}/{run_name}_diagnostics.csv")
if rank == 0:
    with open(f"{outdir}/{run_name}_diagnostics.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "time", "energy", "enstrophy",
            "E_injection", "E_hdiss", "E_drag"])
comm.Barrier()    




##### MAIN LOOP
enstrophy_threshold = 1e10 # stop sim threshold
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        deltaT['g'] = timestep
        solver.step(timestep)
        iterNum['g'] = solver.iteration

        if (solver.iteration-1) % 10 == 0:
            max_psi = np.sqrt(flow.max('psi2'))
            max_u = np.sqrt(flow.max('u2'))
            max_q = np.sqrt(flow.max('q2'))

            # evaluate energetics
            Z  = safe_scalar(1/2*d3.integ(Qtot*Qtot))
            E  = safe_scalar(1/2*d3.integ(u@u + 1/Bu*psi*psi))
            E_inj   = - safe_scalar(d3.integ(psi*forcing))
            E_hdiss = - safe_scalar(nu*d3.integ(psi*d3.lap(a3)))
            E_drag  = safe_scalar(d3.integ(psi*drag))
            residual = E_inj - (E_hdiss + E_drag)

            # write to diagnostics CSV
            with open(f"{outdir}/{run_name}_diagnostics.csv", "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow([
                    solver.sim_time, E, Z,
                    E_inj, E_hdiss, E_drag])
            

            # print status
            logger.info('Iteration=%i, Time=%e, dt=%e, max(psi)=%f, max(u)=%f, max(q)=%f, total enstrophy=%f, total energy=%f' \
                        %(solver.iteration, solver.sim_time,timestep, max_psi, max_u, max_q,Z,E))

            if Z > enstrophy_threshold:
                logger.warning(f"Enstrophy {Z:.3e} exceeded threshold {enstrophy_threshold:.3e} at t={solver.sim_time:.3f}. Exiting loop.")
                break

            elif np.isnan(Z):
                logger.warning(f"Enstrophy is NaN at t={solver.sim_time:.3f}. Exiting loop.")
                break
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()


# plot energetics from CSV file
df = pd.read_csv(f"{outdir}/{run_name}_diagnostics.csv", on_bad_lines='skip')

df = df.sort_values("time")
df = df[df["time"].diff().fillna(1) > 0]

t      = df["time"].to_numpy()
E      = df["energy"].to_numpy()
Z      = df["enstrophy"].to_numpy()
E_inj  = df["E_injection"].to_numpy()
E_hd   = df["E_hdiss"].to_numpy()
E_drag = df["E_drag"].to_numpy()

# running time-average via cumulative trapezoid (smooth, no sawtooth)
def running_trapz_avg(y, t):
    if len(t) < 2:
        return np.full_like(t, np.nan, dtype=float)
    dt = np.diff(t)
    integral = np.concatenate(([0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * dt)))
    return integral / np.maximum(t, np.finfo(float).eps)

E_inj_avg   = running_trapz_avg(E_inj,  t)
E_hdiss_avg = running_trapz_avg(E_hd,   t)
E_drag_avg  = running_trapz_avg(E_drag, t)
residual_avg = E_inj_avg - (E_hdiss_avg + E_drag_avg)

# 3-row subplot: energy, enstrophy, running-avg budget
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(7.0, 9.0), constrained_layout=True)

# Energy
axes[0].plot(t, E, lw=1.8)
axes[0].set_ylabel('Energy')
axes[0].grid(True, alpha=0.35)

# Enstrophy
axes[1].plot(t, Z, lw=1.8)
axes[1].set_ylabel('Enstrophy')
axes[1].grid(True, alpha=0.35)

# Running-averaged budget
axes[2].plot(t, E_inj_avg,   label='Injection')
axes[2].plot(t, E_hdiss_avg, label='HD diss')
axes[2].plot(t, E_drag_avg,  label='Drag diss')
axes[2].plot(t, residual_avg, label='Residual', c='k', lw=2)
axes[2].set_xlabel('Time')
axes[2].set_ylabel('Energy rate (running avg)')
axes[2].grid(True, alpha=0.35)
axes[2].legend()

fig.savefig(f"{outdir}/{run_name}_energetics.png", dpi=200)
