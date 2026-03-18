"""
plot energetics

"""

import csv
import pathlib
import sys
import numpy as np
import dedalus.public as d3
from dedalus.extras import plot_tools
from dedalus.core.domain import *
import os, re, glob
from glob import glob
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # safe on headless clusters
import matplotlib.pyplot as plt
import logging

# plot energetics from CSV file
run_name = sys.argv[1] if len(sys.argv) > 1 else "default"
outdir   = sys.argv[2] if len(sys.argv) > 2 else f"output_{run_name}"
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
diagnostics_csv = os.path.join(outdir, f"{run_name}_diagnostics.csv")
print("Reading file", diagnostics_csv)

# read and clean up time
df = pd.read_csv(diagnostics_csv, on_bad_lines='skip', dtype=str, low_memory=False)
df['time'] = pd.to_numeric(df['time'], errors='coerce')
df = df.dropna(subset=['time'])

# brute-force numeric conversion (coerce any string-like entries)
for col in df.columns:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# now safely extract arrays
t      = df["time"].to_numpy(dtype=float)
E      = df["energy"].to_numpy(dtype=float)
Z      = df["enstrophy"].to_numpy(dtype=float)
E_inj  = df["E_injection"].to_numpy(dtype=float)
E_hd   = df["E_hdiss"].to_numpy(dtype=float)
E_drag = df["E_drag"].to_numpy(dtype=float)

# running time-average via cumulative trapezoid
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


# plot
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(7.0, 6.0), constrained_layout=True)

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

fig.savefig(os.path.join(outdir, f"{run_name}_energetics.png"), dpi=200)
