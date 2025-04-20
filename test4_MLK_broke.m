% Numerically solve for transient conducton problem
clear all; close all; clc;

% Input parameters
L       = 0.5;
H       = 0.5;
dx      = 0.0025;
dy      = 0.0025;
tmax    = 5;
dt      = 0.01;
epsilon = 10^(-4);
alp     = 1.9*10^(-5);
r_x     = alp*dt/dx^2;
r_y     = alp*dt/dy^2;

% This would need to be initialized with WSB, I think
P = 1; % Pressure [atm]
Tign = 300; % Arbitrary ignition temp
Tf = 1000
[rb_wsb, ~]   = project_2_wsb_function(P,Tign); 

% create the x, y meshgrid based on dx, dy
nx    = uint32(L/dx + 1);
ny    = uint32(H/dy + 1);
[X,Y] = meshgrid(linspace(0,L,nx),linspace(0,H,ny));

% Matrix of burning rates - could use this with dt to determine when
% a node burns out. Initialize a burn "start time" for each node, and then
% increment dt until all of that material is burned out. After that, force
% an "off" temperature
% Note this is just initializing, this would be added to the loop
rb = ones(nx,ny);

thickness = 1; % thickness of material [cm]
material = ones(nx,ny).*thickness

% set initial and boundary conditions
T_int = 100;
T     = T_int*ones(ny,nx);

% iteration, march in time
n    = 0; 
nmax = uint32(tmax/dt);

%Logical Matrix
L    = zeros(nx,ny);

% Set corner node to flame temp
T(ny,2) = Tf;
count = 0

while n < nmax

    % Increment time
    n = n + 1;

    T_n = T;

    for j = 2:ny-1
        for i = 2:nx-1
            T(j,i) = T_n(j,i) + r_x*(T_n(j,i+1)-2*T_n(j,i)+T_n(j,i-1))...
                + r_y*(T_n(j+1,i)-2*T_n(j,i)+T_n(j-1,i));
            if L(i,j) == 1 && L_n0(i,j) == 0
               %wsb_test = project_2_wsb_function(P, T(j,i));
               count = count + 1
            end
        end
    end
    
    % for i = 1:nx
    %     j = ny;
    %     T(j,i) = T(j-1, i);
    % end
    % 
    % for i = 1:nx
    %     j = 1;
    %     T(j,i) = T(j+1, i);
    % end
    % 
    % 
    % for j = 1:ny
    %     i = nx;
    %     T(j,i) = T(j, i-1);
    % end
    % 
    % for j = 1:ny
    %     i = 1;
    %     T(j,i) = T(j, i+1);
    % end

    % Search for values exceeding the ignition temperature
    L_n0 = L;
    grid_ignite = (T >= Tign) & (material > 0);
    L = L | grid_ignite;

    % Set any values exceeding ignition temperature to flame temp
    T(L) = Tf;

    % Turn off nodes that have burnt out here, probably?
    burning = L & (material > 0);
    material = material - burning.*(rb*dt);
    material(material < 0) = 0;

    L (material == 0) = 0;


    % if uint16(n/50) == n/50 % refresh the plot every 50 time steps to save time     
    %     contourf(X,Y,T,10);
    %     title(sprintf('Time = %g s',n*dt))
    %     xlabel('x (m)')
    %     ylabel('y (m)')
    %     axis('equal','tight')
    %     pause(0.1)
    % end
    
    if uint16(n/50) == n/50 % refresh the plot every 50 time steps to save time     
        contourf(X,Y,material,0.1);
        title(sprintf('Time = %g s',n*dt))
        xlabel('x (m)')
        ylabel('y (m)')
        axis('equal','tight')
        pause(0.1)
    end

    % check for convergence
    err = max(max(abs((T-T_n))));
    if err <= epsilon
        break
    end
end


