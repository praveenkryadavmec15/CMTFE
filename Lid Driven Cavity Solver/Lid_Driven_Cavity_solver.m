clc
clear 

% Lid Driven Cavity 

rho = 1;
L1 = 1;
L2 = L1;
Lc = L1;

U = 1;
Uc = U; % Characteristic velocity

Re = 100; % Reynolds Number

muu = (rho*Uc*Lc)/Re;

nuu = muu/rho;

e_st = 0.001; % Steady state convergence tolerance
 
imax = 51;
jmax = 51;

Dx = L1/(jmax-1);
Dy = L2/(imax-1);


% Defining Matrix for storing u, v, si and omega value

u = zeros(imax,jmax);
v = zeros(imax,jmax);
si = zeros(imax,jmax);
omega = zeros(imax,jmax);

% Boundary condition

u(:,1) = 0;
u(:,jmax) = 0;
u(imax,:) = 0;
u(1,:) = U;

v(:,1) = 0;
v(:,jmax) = 0;
v(imax,:) = 0;
v(1,:) = 0;

si(:,1) = 0;
si(:,jmax) = 0;
si(imax,:) = 0;
si(1,:) = 0;

omega(:,1) = -2*si(:,2)/Dx^2;
omega(:,jmax) = -2*si(:,jmax-1)/Dx^2;
omega(imax,:) = -2*si(imax-1,:)/Dy^2;
omega(1,:) = -(2*U*Dy + 2*si(2,:))/(Dy^2);

u_old = u;
v_old = v;
si_old = si;
omega_old = omega;

% Iteration count
n = 0;

while true

    n = n + 1;

    %% Calculating u & v velocity %%

    for i = 2:imax-1
        for j = 2:jmax-1

            u(i,j) = (si_old(i-1,j) - si_old(i+1,j))/(2*Dy);

            v(i,j) = -(si_old(i,j+1) - si_old(i,j-1))/(2*Dx);
        end
    end
    

    %% Calculation of streamfunction %%

    for i = 2:imax-1
        for j = 2:jmax-1

            si(i,j) = ((Dy^2) * (si(i,j+1) + si_old(i,j-1)) + (Dx^2) * (si(i+1,j) + si_old(i-1,j)) + (Dx^2)*(Dy^2)* omega_old(i,j))/(2*(Dx^2 + Dy^2));
        end
    end
    

    %% Calculation of omega %% 

    for i = 2:imax-1
        for j = 2:jmax-1

            % coefficient with CD

%             aP = 2*nuu*(1/(Dx^2) + 1/(Dy^2));
%             aE = nuu/Dx^2 - u_old(i,j)/(2*Dx);
%             aW = nuu/Dx^2 + u_old(i,j)/(2*Dx);
%             aN = nuu/Dy^2 - v_old(i,j)/(2*Dy);
%             aS = nuu/Dy^2 + v_old(i,j)/(2*Dy);


            % Coefficient with FOU
            aP = 2*nuu*(1/Dx^2 + 1/Dy^2) + max(u_old(i,j),0)/Dx - min(u_old(i,j),0)/Dx + max(v_old(i,j),0)/Dy - min(v_old(i,j),0)/Dy;
            aE = nuu/Dx^2 - min(u_old(i,j),0)/Dx;
            aW = nuu/Dx^2 + max(u_old(i,j),0)/Dx;
            aN = nuu/Dy^2 - min(v_old(i,j),0)/Dy;
            aS = nuu/Dy^2 + max(v_old(i,j),0)/Dy;

            if u_old(i,j) * Dx/nuu <= 2 && v_old(i,j) * Dy/nuu <= 2
                RF = 1;
            else
                RF = 0.1;
            end

            omega(i,j) = omega_old(i,j) + RF*((aE * omega(i,j+1) + aW * omega_old(i,j-1) + aS * omega(i+1,j) + aN * omega_old(i-1,j))/aP - omega_old(i,j));
                
        end
    end

    % Updating omega BC

    omega(:,1) = -2*si(:,2)/Dx^2;
    omega(:,jmax) = -2*si(:,jmax-1)/Dx^2;
    omega(imax,:) = -2*si(imax-1,:)/Dy^2;
    omega(1,:) = -(2*U*Dy + 2*si(2,:))/(Dy^2);
    

    if max(abs(omega_old - omega)) < e_st & max(abs(si_old - si)) < e_st & max(abs(v_old - v)) < e_st & max(abs(u_old - u)) < e_st
        break
    else
        u_old = u;
        v_old = v;
        si_old = si;
        omega_old = omega;

    end
    %disp(u)
end

disp(n)
disp(u)


%% Plotting of result %%

ucl=u(:,(jmax+1)/2);
vcl=v((imax+1)/2,:);


figure
x=0:Dx:L1;
y=1-(0:Dy:L2);
[X,Y]=meshgrid(x,y);
contourf(X,Y,si,10)
colormap turbo
colorbar
title('Streamlines')
xlabel('X')
ylabel('Y')

figure
x=0:Dx:L1;
y=1-(0:Dy:L2);
[X,Y]=meshgrid(x,y);
contourf(X,Y,u,10)
colormap turbo
colorbar
title('X component of Velocity')
xlabel('X')
ylabel('Y')

figure
x=0:Dx:L1;
y=1-(0:Dy:L2);
[X,Y]=meshgrid(x,y);
contourf(X,Y,v,10)
colormap turbo
colorbar
title('Y component of Velocity')
xlabel('X')
ylabel('Y')

figure
x=0:Dx:L1;
y=1-(0:Dy:L2);
[X,Y]=meshgrid(x,y);
contourf(X,Y,omega,10)
colormap turbo
colorbar
title('Vorticity')
xlabel('X')
ylabel('Y')

figure
y=1-(0:Dy:L2);
y_com = [1;0.9766;0.9688;0.9609;0.9531;0.8516;0.7344;0.6172;0.5;0.4531;0.2813;0.1719;0.1016;0.0703;0.0625;0.0547;0];
u_com=[1;0.84123;0.78871;0.73722;0.68717;0.23151;0.00332;-0.136641;-0.20581;-0.2109;-0.15662;-0.1015;-0.063434;-0.04775;-0.04192;-0.03717;0];
plot(u_com,y_com, 'k *')
xlim([-0.22 1])
hold on
plot(ucl,y, 'g')
ucl31=[1;0.787620129;0.587711723;0.42114596;0.290555896;0.18961776;0.110101319;0.045170674;-0.009748753;-0.057195052;-0.098228157;-0.133016607;-0.161334794;-0.182928957;-0.197740083;-0.206010366;-0.208285895;-0.20535712;-0.198160901;-0.187675667;-0.174825389;-0.160406644;-0.145042783;-0.129160215;-0.112985847;-0.096555039;-0.079727401;-0.062199838;-0.043513792;-0.023049018;0];
ycl31=[1;0.966666667;0.933333333;0.9;0.866666667;0.833333333;0.8;0.766666667;0.733333333;0.7;0.666666667;0.633333333;0.6;0.566666667;0.533333333;0.5;0.466666667;0.433333333;0.4;
0.366666667;0.333333333;0.3;0.266666667;0.233333333;0.2;0.166666667;0.133333333;0.1;0.066666667;0.033333333;0];
plot(ucl31,ycl31,'r')
ucl11=[1;0.400693491;0.087599281;-0.073053304;-0.150595989;-0.168701742;-0.150051734;-0.1161176;-0.079522541;-0.042820845;0];
ycl11=[1;0.9;0.8;0.7;0.6;0.5;0.4;0.3;0.2;0.1;0];
plot(ucl11,ycl11,'b')
legend('Ghia','51x51 grid','31x31 grid','11x11 grid')
title('Centreline velocity')
xlabel('U -->')
ylabel('Y -->')


figure
x=0:Dx:L1;
x_com=[0;0.0625;0.0703;0.0781;0.0938;0.1563;0.2266;0.2344;0.5;0.8047;0.8594;0.9063;0.9453;0.9531;0.9609;0.9688;1];
v_com=[0;0.09233;0.10091;0.1089;0.12317;0.16077;0.17507;0.17527;0.05454;-0.24533;-0.22445;-0.16914;-0.10313;-0.08864;-0.07391;-0.05906;0];
plot(x_com,v_com, 'k *')
hold on
plot(x,vcl, 'g')
vcl31=[0;0.055223461;0.098059308;0.129717495;0.151757886;0.165774834;0.173190206;0.175150795;0.172500499;0.165803711;0.155393962;0.141430839;0.123958774;0.102960071;0.078410224;0.05033074;0.01885264;-0.015716933;-0.052806781;-0.091480194;-0.130327229;-0.167368317;-0.199987923;-0.22496228;-0.238668932;-0.237574833;-0.219059895;-0.182465575;-0.129994098;-0.066823813;0];
xcl31=[0;0.033333333;0.066666667;0.1;0.133333333;0.166666667;0.2;0.233333333;0.266666667;0.3;0.333333333;0.366666667;0.4;0.433333333;0.466666667;0.5;0.533333333;0.566666667;0.6;0.633333333;0.666666667;0.7;0.733333333;0.766666667;0.8;0.833333333;0.866666667;0.9;0.933333333;0.966666667;1];
plot(xcl31,vcl31,'r')
vcl11=[0;0.091991455;0.12281818;0.116607647;0.086399363;0.036138153;-0.031995232;-0.108164468;-0.1616866;-0.136572786;0];
xcl11=[0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1];
plot(xcl11,vcl11,'b')
ylim([-0.26 0.2])
title('Centreline velocity')
legend('Ghia','51x51 grid','31x31 grid','11x11 grid')
xlabel('X -->')
ylabel('V -->')


