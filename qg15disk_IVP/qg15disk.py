"""
Solves the nondimensional 1.5-layer QG equations on a gamma plane

Q = lap(psi) - 1/Bu*psi - 1/Rog*r^2
D/Dt(Q)=0

Burger number Bu = (Ld/L)^2 
Curvature rossby number Rog = 2U/(gamma*L^3) = 2U/(fL)*(ap^2/L^2) = inertial/curvature

Numerics:
    Nondimensionalized with Bu=1 for performance

    want N at least 2*Ldom*N/sqrt(Bu) but a power of 2 for performance
    number of cells per deformation radius N approx 10 if possible, at least 6
    Nr approx Nphi/pi but rounded to a multiple of four

    With parameters below, will run on 32 processes in about 3 hours

To run and plot using e.g. 8 processes:
    
RUN="S2"
DIR="exampleOutput"
mpiexec -n 8 python3 qg15disk.py $RUN $DIR

python3 plot_energetics.py $RUN $DIR

rm -rf frames
mpiexec -n 8 python3 plot_disk_contour.py $RUN $DIR snapshots/*.h5 --output=frames
ffmpeg -y -framerate 24 -i frames/write_%06d.png -c:v libx264 -pix_fmt yuv420p $DIR/${RUN}_output.mp4

"""

import csv
import pathlib
import sys
import numpy as np
from scipy.special import iv  # modified Bessel I_v
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

def zeta_chan(rv, Vm, Rm, b):
    rr = rv / Rm
    return (2*Vm/Rm) * (1 - 0.5*rr**b) * np.exp((1.0/b) * (1 - rr**b))


# configuration
n_vort = 5 # symmetry (n_vort + 1)
gammaplane = 0 # gamma plane (0,1)
BGflow = 1 # background flow (0,1)
# Ld = 1400e3 # estimated NP or SP
# Ld = 1700e3 # NP
Ld = 2000e3 # SP

# parameters from data
# south pole
gamma = 7.7968e-20
Vm_cpc_dim = 86
Rm_cpc_dim = 1050000
Vm_pc_dim  = 88
Rm_pc_dim  = 1300000
b_pc  = 1.35
b = 1.20
r0_cpc = 7010000 # from the precession dataset

# # north pole
# gamma = 7.7968e-20
# Vm_cpc_dim = 44
# Rm_cpc_dim = 1100000
# Vm_pc_dim  = 75
# Rm_pc_dim  = 900000
# b_pc  = 1.15
# b = 1.10
# r0_cpc = 8000000 # from the precession dataset

# # modify for stability
# Rm_pc_dim = 1100000
# Vm_pc_dim  = 100
# b_pc  = 1.0
# b = 2.0

# nondimensionalization
U = 100 
L = Ld
Bu = (Ld/L)**2
Rog = 2*U/(gamma*L**3)
Vm = Vm_cpc_dim/U
Rm = Rm_cpc_dim/L
Vm_pc = Vm_pc_dim/U
Rm_pc = Rm_pc_dim/L
Rprec = r0_cpc/L
Rdom = 3*Rprec

# Simulation parameters
sim_time_years = 10
stop_sim_time = sim_time_years*3600*24*365*U/L
nu = 1e-10
Nphi = 512
Nr = int(np.round(Nphi/np.pi/4)*4)
dealias = 3/2
write_speed = stop_sim_time/1000
timestepper = d3.RK222
max_timestep = 1e-1
dtype = np.float64

# forcing
ell = 8*np.pi*Rdom/Nphi # forcing scale in physical space
Qpert = 0.01 # random noise perturbation strength

# Bases
coords = d3.PolarCoordinates('phi', 'r')
dist = d3.Distributor(coords, dtype=dtype) # change mesh based on number of processes
disk = d3.DiskBasis(coords, shape=(Nphi, Nr), radius=Rdom, dealias=dealias, dtype=dtype)
edge = disk.edge
domain = Domain(dist, bases=[disk])

# Fields
t = dist.Field()
deltaT = dist.Field()
iterNum = dist.Field()
eta = dist.Field(name='eta', bases=disk)
noise = dist.Field(name='noise', bases=disk)
Ubg = dist.VectorField(coords, bases=disk, name='Ubg')
Qbg = dist.Field(bases=disk, name='Qbg')

psi = dist.Field(name='psi', bases=disk)
q0 = dist.Field(name='q0',bases=disk)
qf = dist.Field(name='qf', bases=disk)
q_vort = dist.Field(name='q_vort',bases=disk)
zeta_vort = dist.Field(name='zeta_vort',bases=disk)
psi_vort = dist.Field(name='psi_vort', bases=disk)

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
dr = lambda A: d3.radial(d3.grad(A)(r=Rdom))
ephi = dist.VectorField(coords, bases=disk)
ephi['g'][0], ephi['g'][1] = 1, 0
u = d3.Skew(d3.Gradient(psi))
HD = nu*d3.lap(a3) # hyperdiffusion (4th order)

# gamma plane case
if gammaplane==1:
    eta['g'] = -1/Rog*r**2

# background flow
if BGflow==1:
    eta['g'] = 0
    Rmix_hat = Rdom
    Rd_hat = np.sqrt(Bu)
    Qbg['g'] = -1/Rog*(1/2*Rmix_hat**2)
    Ubg['g'][0] = - 2/Rog*Bu * (r - (iv(1, r / Rd_hat) / iv(1, Rmix_hat / Rd_hat)) * Rmix_hat)
    Ubg['g'][1] = 0 


## BUILD SOLVERS
# build helmholtz filter for noise
problem = d3.LBVP([qf, a1, a2, a3, tau_q, tau_a1, tau_a2, tau_a3],namespace=locals())
problem.add_equation("a1 - lap(qf)  + lift(tau_a1) = 0")
problem.add_equation("a2 - lap(a1)  + lift(tau_a2) = 0")
problem.add_equation("a3 - lap(a2)  + lift(tau_a3) = 0")
problem.add_equation("qf + ell**8*lap(a3) + lift(tau_q) = q0")
problem.add_equation("qf(r=Rdom) = 0")
problem.add_equation("a1(r=Rdom) = 0")
problem.add_equation("a2(r=Rdom) = 0")
problem.add_equation("a3(r=Rdom) = 0")
forcing_filter = problem.build_solver()

# calculate psi from zeta
psi_invert_lbvp = d3.LBVP([psi_vort, tau_psi], namespace=locals())
psi_invert_lbvp.add_equation("lap(psi_vort) + lift(tau_psi) = zeta_vort")
psi_invert_lbvp.add_equation("psi_vort(r=Rdom) = 0")
psi_invert = psi_invert_lbvp.build_solver()

# Main problem
problem = d3.IVP([q, psi, a1, a2, a3, tau_q, tau_a1, tau_a2, tau_a3, tau_psi],time=t, namespace=locals())
problem.add_equation("dt(q) + HD + lift(tau_q) = - (u + Ubg)@grad(q + eta)")
problem.add_equation("a1 - lap(q)  + lift(tau_a1) = 0")
problem.add_equation("a2 - lap(a1) + lift(tau_a2) = 0")
problem.add_equation("a3 - lap(a2) + lift(tau_a3) = 0")
problem.add_equation("q - lap(psi) + 1/Bu*psi + lift(tau_psi) = 0")

problem.add_equation("psi(r=Rdom) = 0") # material surface
problem.add_equation("q(r=Rdom)  = 0") # free slip
problem.add_equation("a1(r=Rdom)  = 0")
problem.add_equation("a2(r=Rdom)  = 0")
problem.add_equation("a3(r=Rdom)  = 0")

# Main solver setup
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time
file_handler_mode = 'overwrite'
initial_timestep = max_timestep/10


##### INITIAL CONDITIONS
q.change_scales(dealias)
q0.change_scales(dealias)
qf.change_scales(dealias)

# random noise
q0.fill_random('g', seed=42, distribution='normal', scale=1.0)
forcing_filter.solve()
Qrms = np.sqrt(1/(np.pi*Rdom**2) * safe_scalar(d3.integ(qf*qf)))
qf['g'] *= Qpert / Qrms
q['g'] += qf['g']

# Add circumpolar vortices
for i in range(n_vort):
    phi0 = i * 2*np.pi / n_vort
    r0 = Rprec
    rv = np.sqrt(r**2 + r0**2 - 2*r*r0*np.cos(phi - phi0))
    zeta_vort['g'] += zeta_chan(rv, Vm, Rm, b)

# Add central vortex
rv = r
zeta_vort['g'] += zeta_chan(rv, Vm_pc, Rm_pc, b_pc)

# calculate q and add to initial condition
psi_invert.solve()
zeta_vort.change_scales(dealias)
psi_vort.change_scales(dealias)
q['g'] += zeta_vort['g'] - 1/Bu*psi_vort['g']


# Analysis
Qtot = Qbg + q + eta # total PV
zeta = d3.lap(psi)
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=write_speed, max_writes=10, mode=file_handler_mode)
snapshots.add_task(psi, name='psi')
snapshots.add_task(zeta, name='zeta')
snapshots.add_task(q, name='q')
snapshots.add_task(Qtot, name='Qtot')
snapshots.add_task(u, name='u')
snapshots.add_task(u + Ubg, name='Utot')
snapshots.add_task(np.sqrt(u@u), name='speed')

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
    "Vm_cpc": Vm_cpc_dim,
    "Rm_cpc": Rm_cpc_dim,
    "Vm_pc": Vm_pc_dim,
    "Rm_pc": Rm_pc_dim,
    "b_pc": b_pc,
    "b_cpc": b,
    "Rprec": r0_cpc,
    "Ld": Ld,
    "Runtime years": sim_time_years,
    "Qpert": Qpert,
    "nu": nu,
    "Nphi": Nphi,
    "Nr": Nr,
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
            "time", "energy", "enstrophy", "E_hdiss"])
comm.Barrier()    


##### MAIN LOOP
enstrophy_threshold = 1e10 # stop sim threshold
print('Stop_sim_time = ',stop_sim_time)
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
            E_hdiss = - safe_scalar(nu*d3.integ(psi*d3.lap(d3.lap(d3.lap(d3.lap(q))))))

            # write to diagnostics CSV
            with open(f"{outdir}/{run_name}_diagnostics.csv", "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow([
                    solver.sim_time, E, Z,E_hdiss])
            

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