# VortexCrystal
Simple models for vortex interactions

1) The matlab code contourDynamics.m simulates vortex patch interaction with two examples of ring vortices. The first example is stable, and vortices propagate under mutual advection without merging. In the second example vortices are unstable and merge via deformation. RK2 timestepping is used to better conserve vortex patch area. The contour dynamics code does NOT do any surgery or remeshing, and instead has a simple area cutoff for stability. 

