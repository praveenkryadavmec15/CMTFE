**2D Incompressible Flow steady state solver for a Lid Driven Cavity and a T-channel are provided here**  
These codes were written as part of the Term project for Computational Methods in Thermal and Fluids Engineering course at IIT Bombay.

**Domain Discretization: Grid Generation**  
Uniform grid was used for dividing the geometry into smaller grids.

**Finite Difference Method Based Algebraic Formulation**  
Finite Difference Method (FDM) was used for the discretization of the Governing Differential Equations (stream function and vorticity equation) to Linear Algebraic Equations.   
1. Finite Difference Equation for Stream Function
2. Finite Difference Equation for u and v velocity
3. Finite Difference Equation for Vorticity Transport Equation

**Solution Algorithm**  
1. Input Parameters
   1. $i_{max}$ and $j_{max}$ are specified. Matrix for $\phi$ ( = si, omega, u & v) and $\phi_{old}$ having $i_{max}$ x $j_{max}$ elements is specified.
   2. $U_o$ , L & $\epsilon$ are specified.
   3. Reynolds number specified and $\Delta X$ and $\Delta Y$ are calculated.
2. Initial value for  $\phi$ ( = si, omega, u & v) is assumed.
3. Boundary condition for si, u and v are applied.
4. $\phi_{old}$(i,j) = $\phi_{old}$(i,j)  $\phi$ ( = si, omega, u & v) for all grid points.
5. Velocities are computed at the interior nodes using stream function values of previous iteration.
6. Stream function values at interior nodes are computed using vorticity values of previous iteration.
7. Vorticity values at the domain boundaries are computed.
8. Vorticity values at the interior nodes are computed.
9. For strongly nonlinear problems (high Reynolds number flows), under relaxation may be employed for si and omega.
10. For each node, the difference in the values of u, v, si and omega between two consecutive iterations are determined. Maximum residual value is determined among them. If the residue > $\epsilon$ one goes back to step 4.

**CFD Analysis**  
u & v velocity contour, streamlines and vorticity contours have been plotted. 
