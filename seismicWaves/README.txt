This folder contains two codes for numerically computing a seismic waves in rectangular domains. The two codes have different boundary conditions for testing the strength of the solver algorithm. When computing, I have used first order finite elements for the spatial domains, and finite differences for the time domains. This ensures second order errors. I have then run convergence tests with P and S waves to validate the solutions. The principal differnces in the solvers are the boundary conditions at the top of the rectangles. 

seismic-wave-dirichlet: Dirichlet BC-s on all sides
seismic-wave-stress   : Free surface BC on the top

NOTE: The solvers are still under development, and will be generalized to classes in the near future.
