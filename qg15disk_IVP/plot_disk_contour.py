"""
Plot disk outputs.

Usage:
    plot_disk.py <run_name> <outdir> <files>... [--output=<dir>]

Options:
    --output=<dir>  Output directory [default: ./frames]
"""

import pathlib
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main(filename, start, count, output):
    """Save plot of specified tasks for given range of analysis writes."""

    task = "Qtot"
    clim_set = (-2, 2)
    #clim_set = (0, 0) # flexible bounds
    cmap = "RdBu_r"
    dpi = 200
    nlevels = 21

    def savename_func(write):
        return f"write_{write:06}.png"

    def title_func(sim_time):
        return f"t = {sim_time:.3f}"

    with h5py.File(filename, mode="r") as file:
        dset = file["tasks"][task]

        # Assumes dataset shape is (time, phi, r)
        phi_name = dset.dims[1][0].name
        r_name = dset.dims[2][0].name

        phi = file[phi_name][...]
        r = file[r_name][...]

        for index in range(start, start + count):

            data = dset[index, :, :]

            # 1. Close azimuthal seam
            phi_plot = np.concatenate([phi, [phi[0] + 2*np.pi]])
            data_plot = np.vstack([data, data[0:1, :]])

            # 2. Add r = 0 to remove central hole
            r_plot = np.concatenate([[0.0], r])

            # Use azimuthal mean at inner radius
            center_value = np.mean(data_plot[:, 0])
            center_col = np.full((data_plot.shape[0], 1), center_value)

            data_plot = np.hstack([center_col, data_plot])

            # 3. Build mesh using modified coordinates
            phi_mesh, r_mesh = np.meshgrid(phi_plot, r_plot, indexing="ij")
            x_mesh = r_mesh * np.cos(phi_mesh)
            y_mesh = r_mesh * np.sin(phi_mesh)

            if clim_set == (0, 0):
                lim = np.nanmax(np.abs(data_plot))
                levels = np.linspace(-lim, lim, nlevels)
            else:
                levels = np.linspace(clim_set[0], clim_set[1], nlevels)

            fig, ax = plt.subplots(figsize=(5, 6), constrained_layout=False)
            fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.8)

            cf = ax.contourf(
                x_mesh,
                y_mesh,
                data_plot,
                levels=levels,
                cmap=cmap,
                extend="both",
            )

            # # Solid contour lines
            # # line_levels = np.linspace(levels[0], levels[-1], 21)
            # line_levels = levels
            # ax.contour(
            #     x_mesh,
            #     y_mesh,
            #     data_plot,
            #     levels=line_levels,
            #     colors="k",
            #     linewidths=0.4,
            #     linestyles="solid",   # force solid lines
            #     alpha=0.6,
            # )

            domain_radius = 0.5*float(np.max(r_plot))
            ax.set_aspect("equal")
            ax.set_xlim(-domain_radius, domain_radius)
            ax.set_ylim(-domain_radius, domain_radius)
            ax.axis("off")

            # --- Colorbar at top ---
            cbar = fig.colorbar(
                cf,
                ax=ax,
                orientation="horizontal",
                location="top",   # <-- key change
                pad=0.05,         # spacing from axes
                fraction=0.05,
            )
            cbar.outline.set_visible(False)
            cbar.set_label(task)

            # --- Title above colorbar ---
            sim_time = file["scales/sim_time"][index]
            fig.suptitle(
                title_func(sim_time),
                x=0.5,
                y=0.95,   # keep safely above colorbar
                ha="center",
            )

            write_number = file["scales/write_number"][index]
            savepath = output.joinpath(savename_func(write_number))

            # Important: do not use bbox_inches="tight"
            fig.savefig(savepath, dpi=dpi)
            plt.close(fig)


if __name__ == "__main__":
    from docopt import docopt
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    run_name = args["<run_name>"]
    outdir = args["<outdir>"]

    output_path = pathlib.Path(args["--output"]).absolute()

    with Sync() as sync:
        if sync.comm.rank == 0 and not output_path.exists():
            output_path.mkdir()

    post.visit_writes(args["<files>"], main, output=output_path)