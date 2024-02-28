clc
clear 

% Flow in a T shaped plane channel with one inlet and two outlet 

rho = 1; 
L1 = 15;
L2 = 5;
H = 1; 
Lc = H;

V = 1; % Flow velocity at inlet
Uc = V; % Characteristic Velocity at inlet

Re = input("Enter Reynolds Number:"); % Reynolds Number

muu = (rho*Uc*Lc)/Re;

nuu = muu/rho;

e_st = 0.001; % Steady state convergence tolerance

imax = 151;
jmax = 151;

Dx = L1/(imax-1); % x-dir grid size
Dy = (L2 + H)/(jmax-1); % y-dir grid size

% Defining matrix for u,v,si and omega

u = zeros(imax,jmax);
v = zeros(imax,jmax);
si = zeros(imax,jmax);
omega = zeros(imax,jmax);

% Boundary condition

u(1,127:150) = u(2,127:150); % Left Outlet
u(imax,127:150) = u(imax-1,127:150); % Right Outlet

u(72:80,1) = 0; % Bottom Inlet 

u(:,jmax) = 0; % Top Wall

u(1:71,126) = 0; % Left Horizontal Wall
u(71,1:126) = 0; % Left Vertical Wall

u(81:151,126) = 0; % Right Horizontal Wall
u(81,1:126) = 0; % Right Vertical Wall


v(1,127:150) = 0; % Left Outlet
v(imax,127:150) = 0; % Right Outlet

v(72:80,1) = V; % Bottom Inlet

v(:,jmax) = 0; % Top Wall

v(1:71,126) = 0; % Left Horizontal Wall
v(71,1:126) = 0; % Left Vertical Wall

v(81:151,126) = 0; % Right Horizontal Wall
v(81,1:126) = 0; % Right Vertical Wall


si(1,127:150) = si(2,127:150); % Left Outlet
si(imax,127:150) = si(imax-1,127:150); % Right Outlet

si(72:80,1) = si(72:80,2); % Bottom Inlet

si(:,jmax) = V*H/2; % Top wall

si(1:71,126) = V*H; % Left Horizontal Wall
si(71,1:126) = V*H; % Left Vertical Wall

si(81:151,126) = 0; % Right Horizontal Wall
si(81,1:126) = 0; % Right Vertical Wall


omega(1,127:150) = -(si(1,128:151) -2*si(1,127:150) + si(1,126:149))/Dy^2; % Left Outlet
omega(imax,127:150) = -(si(imax,128:151) -2*si(imax,127:150) + si(imax,126:149))/Dy^2; % Right Outlet 

omega(72:80,1) = -(si(73:81,1) -2*si(72:80,1) + si(71:79,1))/Dx^2; % Bottom Inlet

omega(:,jmax) = (2*si(:,jmax) -2*si(:,jmax-1) -2*u(:,jmax)*Dy)/(Dy^2); % Top wall

omega(1:71,126) = (2*si(1:71,126) -2*si(1:71,127) +2*u(1:71,126)*Dy)/(Dy^2); % Left Horizontal Wall
omega(71,1:126) = (2*si(71,1:126) -2*si(72,1:126) -2*v(71,1:126)*Dx)/(Dx^2); % Left Vertical Wall

omega(81:151,126) = (2*si(81:151,126) -2*si(81:151,127) + 2*u(81:151,126)*Dy)/(Dy^2); % Right Horizontal Wall
omega(81,1:126) = (2*si(81,1:126) -2*si(80,1:126) +2*v(81,1:126)*Dx)/(Dx^2); % Right Vertical Wall


u_old = u;
v_old = v;
si_old = si;
omega_old = omega;

% Iteration Count
n = 0;

while true

    n = n + 1;
 
    %% Boundary Conditions %%

    % u-velocity BC

    u(1,127:150) = u_old(2,127:150); % Left Outlet
    u(imax,127:150) = u_old(imax-1,127:150); % Right Outlet
    
    u(72:80,1) = 0; % Bottom Inlet 
    
    u(:,jmax) = 0; % Top Wall
    
    u(1:71,126) = 0; % Left Horizontal Wall
    u(71,1:126) = 0; % Left Vertical Wall
    
    u(81:151,126) = 0; % Right Horizontal Wall
    u(81,1:126) = 0; % Right Vertical Wall

    % v-velocity BC

    v(1,127:150) = 0; % Left Outlet
    v(imax,127:150) = 0; % Right Outlet
    
    v(72:80,1) = V; % Bottom Inlet
    
    v(:,jmax) = 0; % Top Wall
    
    v(1:71,126) = 0; % Left Horizontal Wall
    v(71,1:126) = 0; % Left Vertical Wall
    
    v(81:151,126) = 0; % Right Horizontal Wall
    v(81,1:126) = 0; % Right Vertical Wall


    % Si BC

    si(1,127:150) = si_old(2,127:150); % Left Outlet
    si(imax,127:150) = si_old(imax-1,127:150); % Right Outlet
    
    si(72:80,1) = si_old(72:80,2); % Bottom Inlet
    
    si(:,jmax) = V*H/2; % Top wall
    
    si(1:71,126) = V*H; % Left Horizontal Wall
    si(71,1:126) = V*H; % Left Vertical Wall
    
    si(81:151,126) = 0; % Right Horizontal Wall
    si(81,1:126) = 0; % Right Vertical Wall


    % Omega BC 

    omega(1,127:150) = -(si_old(1,128:151) -2*si_old(1,127:150) + si_old(1,126:149))/Dy^2; % Left Outlet
    omega(imax,127:150) = -(si_old(imax,128:151) -2*si_old(imax,127:150) + si_old(imax,126:149))/Dy^2; % Right Outlet 
    
    omega(72:80,1) = -(si_old(73:81,1) -2*si_old(72:80,1) + si_old(71:79,1))/Dx^2; % Bottom Inlet
    
    omega(:,jmax) = (2*si_old(:,jmax) -2*si_old(:,jmax-1) -2*u_old(:,jmax)*Dy)/(Dy^2); % Top wall
    
    omega(1:71,126) = (2*si_old(1:71,126) -2*si_old(1:71,127) +2*u_old(1:71,126)*Dy)/(Dy^2); % Left Horizontal Wall
    omega(71,1:126) = (2*si_old(71,1:126) -2*si_old(72,1:126) -2*v_old(71,1:126)*Dx)/(Dx^2); % Left Vertical Wall
    
    omega(81:151,126) = (2*si_old(81:151,126) -2*si_old(81:151,127) + 2*u_old(81:151,126)*Dy)/(Dy^2); % Right Horizontal Wall
    omega(81,1:126) = (2*si_old(81,1:126) - 2*si_old(80,1:126) + 2*v_old(81,1:126)*Dx)/(Dx^2); % Right Vertical Wall

    %% Calculation of u & v velocity %%

    for i = 2:150
        for j = 127:150
            u(i,j) = (si_old(i,j+1) - si_old(i,j-1))/(2*Dy);
            v(i,j) = -(si_old(i+1,j) - si_old(i-1,j))/(2*Dx);
        end
    end

    for i = 72:80
        for j = 2:126
            u(i,j) = (si_old(i,j+1) - si_old(i,j-1))/(2*Dy);
            v(i,j) = -(si_old(i+1,j) - si_old(i-1,j))/(2*Dx);
        end            
    end

    %% Calculation of streamfunction %%

    Aw = 1/Dx^2;
    Ae = Aw;
    An = 1/Dy^2;
    As = An;
    Ap = Aw + Ae + An + As;

    for i = 2:150
        for j = 127:150
            si(i,j) = (omega_old(i,j) + Aw*si_old(i-1,j) + Ae*si_old(i+1,j) + As*si_old(i,j-1) + An*si_old(i,j+1))/Ap ;
        end
    end

    for i = 72:80
        for j = 2:126
            si(i,j) = (omega_old(i,j) + Aw*si_old(i-1,j) + Ae*si_old(i+1,j) + As*si_old(i,j-1) + An*si_old(i,j+1))/Ap ;
        end
    end


    %% Calculation of vorticity %%

    for i = 2:150
        for j = 127:150

            RF = 0.1  ;

            aE = nuu/Dx^2 - min(u_old(i,j),0)/Dx;
            aW = nuu/Dx^2 + max(u_old(i,j),0)/Dx;
            aN = nuu/Dy^2 - min(v_old(i,j),0)/Dy;
            aS = nuu/Dy^2 + max(v_old(i,j),0)/Dy;

            aP = aE + aW + aN + aS;

            omega(i,j) = omega_old(i,j) + RF*((aE*omega_old(i+1,j) + aW*omega_old(i-1,j) + aS*omega_old(i,j-1) + aN*omega_old(i,j+1))/aP - omega_old(i,j));
        end    
    end
    
    for i = 72:80
        for j = 2:124

            RF = 0.1;

            aE = nuu/Dx^2 - min(u_old(i,j),0)/Dx;
            aW = nuu/Dx^2 + max(u_old(i,j),0)/Dx;
            aN = nuu/Dy^2 - min(v_old(i,j),0)/Dy;
            aS = nuu/Dy^2 + max(v_old(i,j),0)/Dy;
            aP = aE + aW + aN + aS;
             
    
            omega(i,j) = omega_old(i,j) + RF*((aE*omega_old(i+1,j) + aW*omega_old(i-1,j) + aS*omega_old(i,j-1) + aN*omega_old(i,j+1))/aP - omega_old(i,j));
        end    
    end


    %% Checking for Convergence %%

    if max(abs(omega_old - omega)) < e_st & max(abs(si_old - si)) < e_st & max(abs(v_old - v)) < e_st & max(abs(u_old - u)) < e_st
        break
    else
        u_old = u;
        v_old = v;
        si_old = si;
        omega_old = omega;

    end
    %disp(si)
end

disp(n)
disp(u)

%% Tagging the point outside the geometrical domain as NaN %%

for i = 1:70
    for j = 1:125
         si (i,j)= NaN;
         omega (i,j)= NaN;
         u (i,j)= NaN;
         v (i,j)= NaN;
    end
end
for i = 82:151
    for j = 1:125
         si (i,j)= NaN;
         omega (i,j)= NaN;
         u (i,j)= NaN;
         v (i,j)= NaN;
    end
end
    
%% Plotting of results %%

figure

x = 0:Dx:L1;
y = 0:Dy:L2+H;

[X,Y] = meshgrid((6/15)*x,(15/6)*y);
contourf(Y,X,omega,10)
colormap jet
colorbar

title("Vorticity")
xlabel("X")
ylabel("Y")

figure

x = 0:Dx:L1;
y = (0:Dy:L2+H);

[X,Y] = meshgrid((6/15)*x,(15/6)*y);
contourf(Y,X,si,10)
colormap jet
colorbar

title("Streamlines")
xlabel("X")
ylabel("Y")

figure

x = 0:Dx:L1;
y = (0:Dy:L2+H);

[X,Y] = meshgrid((6/15)*x,(15/6)*y);
contourf(Y,X,u,10)
colormap jet
colorbar

title("x-comp of velocity (u)")
xlabel("X")
ylabel("Y")

figure

x = 0:Dx:L1;
y = (0:Dy:L2+H);

[X,Y] = meshgrid((6/15)*x,(15/6)*y);
contourf(Y,X,v,10)
colormap jet
colorbar

title("y-comp of velocity (v)")
xlabel("X")
ylabel("Y")
