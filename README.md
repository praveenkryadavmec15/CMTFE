**2D Incompressible Flow steady state solver for a Lid Driven Cavity and a T-channel are provided here**

These codes were written as part of the Term project for Computational Methods in Thermal and Fluids Engineering course at IIT Bombay.

**Domain Discretization: Grid Generation**

Uniform grid was used for the dividing the geometry into smaller grids.

**Finite Difference Method Based Algebraic Formulation**

Finite Difference Method (FDM) was used for the discretization of the Governing Differential Equations (stream function and vorticity equation) to Linear Algebraic Equations. 
    1. Finite Difference Equation for Stream Function
    2. Finite Difference Equation for u and V velocity
    3. Finite Difference Equation for Vorticity Transport Equation


The governing differential equation was solved on each grid point obatiaon after geometry discretization to obtain the flow properties 
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
