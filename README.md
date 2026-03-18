# VortexCrystal
Simple models for vortex interactions motivated by Juno observations of the polar vortex crystals on Jupiter, see https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022JE007241 

1) To study the motion of the polar cyclones, contourDynamics.m simulates vortex patch interaction with two examples of ring vortices. The first example is stable, and vortices propagate under mutual advection without merging. In the second example vortices are unstable and merge via deformation. RK2 timestepping is used to better conserve vortex patch area. The contour dynamics code does NOT do any surgery or remeshing, and instead has a simple area cutoff for stability. 

2) Unifying the formation and motion of the cyclones is an area of active study. The Dedalus model qg15diskFD.py simulates the 1.5 layer quasi-geostrophic equations with Ornstein-Uhlenbeck forcing and quadratic drag dissipation, mimicing convective forcing at the polar regions. In the example provided, energy moves up scale to a dipole structure and no vortex crystal emerges, but the numerical setup seems promising for further exploration of the parameter space. Plotting and diagnostic scripts are included.
