A 2D incompressible steady state solver for a Lid Driven Cavity Flow has been written in MATLAB. 

The Geometric and Thermophysical properties of the fluid flow were specified.  
The Geometry is discretized using uniform grid spacing in the x and y direction.
The Boundary conditions and steady state convergence tolerance were also specified.
The flow variable u, v, si (streamfunction), omega (Vorticity) were initialized at each grid point.
The Boundary grid points were assigned values based on the BCs (Neuman or Dirichlet BCs)
Values of all the flow variables were specified or calculated at the boundary based on the boundary condition (Neuman or Dirichlet).
For the inner grid points u and v velocities are obtained from si values.
si values are further updated using old si and omega values.
omega values are updated using old values of omega, u and v velocities.
The omega value for boundary grid points are also updated using the new si values.
If the updated value of flow variables u, v, si and omega doesn't differ with the old values by more than the steady state convergence tolerace the code is stopped
Else the updated values are used as the old values for the next iteration.
The plots for all the flow variables have been plotted.
