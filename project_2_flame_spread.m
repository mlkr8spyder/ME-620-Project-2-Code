% Numerically solve for transient conducton problem
clear all; close all; clc;

% Input parameters
L       = 0.3;
H       = 0.05;
dx      = 0.001;
dy      = 0.0005;
tmax    = 5;
dt      = 0.0005;
epsilon = 10^(-4);
k = 0.2;
cp = 1400;
rho = 1800;
% alp = k / (cp * rho);
alp     = 1.9*10^(-5);
r_x     = alp*dt/dx^2;
r_y     = alp*dt/dy^2;

P = 5; % Pressure [atm]
Tign = 500; % Arbitrary ignition temp
Tf = 2000
[rb_wsb, Tf_test, T_flame]   = project_2_wsb_function(P,Tign); 

% Tf = Tf;
% % Tf = T_flame

% Tf = Tf;
Tf = T_flame

% create the x, y meshgrid based on dx, dy
nx    = uint32(L/dx + 1);
ny    = uint32(H/dy + 1);
[X,Y] = meshgrid(linspace(0,L,nx),linspace(0,H,ny));

thickness = 0.1; % thickness of material [cm]
material = ones(ny,nx).*thickness;

% set initial and boundary conditions
T_int = 293;
T = T_int*ones(ny,nx);

% iteration, march in time
n    = 0; 
nmax = uint32(tmax/dt);

%Logical Matrix
L = zeros(ny,nx);

% Set corner node to flame temp
% T(floor(ny/2),1) = Tf;
T((ny*0.4):(ny*0.6),2) = Tf;

count = 0;

grid_ignite = (T >= Tign) & (material > 0);
L = L | grid_ignite;
burning = L & (material > 0);
L_n0 = L;

test_time_history = zeros(3,nmax);

% while any(burning(:)) && n < nmax
while n < nmax
    
    % Increment time
    n = n + 1;
    
    T_n = T;

    for j = 2:ny-1
        for i = 2:nx-1
            T(j,i) = T_n(j,i) + r_x*(T_n(j,i+1)-2*T_n(j,i)+T_n(j,i-1))...
                + r_y*(T_n(j+1,i)-2*T_n(j,i)+T_n(j-1,i));
            if L(j,i) == 1 && L_n0(j,i) == 0
               %wsb_test = project_2_wsb_function(P, T(j,i));
               count = count + 1;
            end
        end
    end

    test_j = 25;
    test_i = 10;

    test_time_history(:,n) = [T(test_j,test_i); burning(test_j,test_i); material(test_j,test_i)];
    
    % Search for values exceeding the ignition temperature
    L_n0 = L;
    grid_ignite = (T >= Tign) & (material > 0);
    L = L | grid_ignite;

    % Set any values exceeding ignition temperature to flame temp
    T(L) = Tf;

    % Turn off nodes that have burnt out here, probably?
    burning = L & (material > 0);
    material = material - burning.*(rb_wsb*dt);
    material(material < 0) = 0;

    L (material == 0) = 0;


    plot_refresh = 50;
    if uint16(n/plot_refresh) == n/plot_refresh % refresh the plot every 50 time steps to save time     
        clf
        
        subplot(1,3,1)
        contourf(X, Y, T,  10, 'LineColor','none')
        title(sprintf('T @ %g s', n*dt))
        xlabel('x (m)')
        ylabel('y (m)')
        colorbar
        axis equal tight
        
        subplot(1,3,2)
        contourf(X, Y, material, 10, 'LineColor','none')
        title(sprintf('Thickness @ %g s', n*dt))
        xlabel('x (m)')
        ylabel('y (m)')
        axis equal tight
        
        subplot(1,3,3)
        contourf(X, Y, L, 2, 'LineColor','none')
        title(sprintf('Burning @ %g s', n*dt))
        xlabel('x (m)')
        ylabel('y (m)')
        axis equal tight
        
        pause(0.1)
    end
    

    % check for convergence
    err = max(max(abs((T-T_n))));
    if err <= epsilon
        break
    end
end

test_time_history = test_time_history(:,1:n);

%% Time History Plotting
figure
plot(1:n,test_time_history(1,:))
ylim([0 2200])
hold on
grid on
yyaxis right
plot(1:n,test_time_history(3,:))

