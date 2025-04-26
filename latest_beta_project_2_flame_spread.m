% Numerically solve for transient conducton problem
clear all; close all; clc;

% Input parameters
grid_length = 10;    % cm  
grid_height = 5;     % cm
dx = 1.0;            % cm  (1 cm spacing)
dy = 0.5;            % cm  (0.5 cm spacing)
time_max = 10;
dt = 0.005;
tolerance = 10^(-4);
k = 0.2;
cp = 1400;
rho = 1800;
alpha   = (k/(cp*rho)) * 1e4;
r_x     = alpha*dt/dx^2;
r_y     = alpha*dt/dy^2;

k_gas   = 0.07;       % W/(m·K), e.g. air
cp_gas  = 1005;        % J/(kg·K)
rho_gas = 1.2;         % kg/m^3
alpha_g = (k_gas/(cp_gas*rho_gas)) * 1e4;
r_x_g   = alpha_g*dt/dx^2;
r_y_g   = alpha_g*dt/dy^2;

P = 1; % Pressure [atm]
temp_ign = 400; % Arbitrary ignition temp
temp_flame = 2000
[rb_wsb, Tf_test, temp_flame_wsb]   = project_2_wsb_function(P,temp_ign); 

% create the x, y meshgrid based on dx, dy
nx = uint32(grid_length/dx + 1);
ny = uint32(grid_height/dy + 1);
[X,Y] = meshgrid(linspace(0,grid_length,nx),linspace(0,grid_height,ny));

thickness = 0.5; % thickness of material [cm]
material = ones(ny,nx).*thickness;

% set initial and boundary conditions
temp_int = 293;
temp = temp_int*ones(ny,nx);

% iteration, march in time
n = 0; 
nmax = uint32(time_max/dt);

%Logical Matrix
node_status = zeros(ny,nx);

% Set corner node to flame temp
% T(floor(ny/2),1) = Tf;
temp((ny*0.4):(ny*0.6),2) = temp_flame;
material(1,  :) = 0;    % bottom row
material(end,:) = 0;    % top row
material(:,  1) = 0;    % left  column
material(:,end) = 0;    % right column

%–– initialize local burn‐rate & flame‐temp maps
rb_map = zeros(ny,nx);
Tf_map = temp_int * ones(ny,nx);

%–– seed the spark region with its WSB values
js = round(ny*0.4):round(ny*0.6);
is = 2;
[rb0, ~, Tf0] = project_2_wsb_function(P, temp_int);
rb_map(js,is)   = rb0;   % cm/s
Tf_map(js,is)   = Tf0;   % K

count = 0;

grid_ignite = (temp >= temp_ign) & (material > 0);
node_status = node_status | grid_ignite;
burning = node_status & (material > 0);
node_status_n0 = node_status;

test_time_history = zeros(3,nmax);
time_history_j = round(ny/2);
time_history_i = round(nx/2);

frontPos  = zeros(1,nmax);    % store front location (cm)
frontRate = zeros(1,nmax);    % store instantaneous spread rate (cm/s)


[~, Tsd_amb, ~] = project_2_wsb_function(P, temp_int);
temp_ign = Tsd_amb;    % real ignition temp [K]

%–– build X-implicit matrix for cols 2…nx−1
ex = ones(nx-2,1);
Ax = spdiags([ -r_x*ex, (1+2*r_x)*ex, -r_x*ex ], -1:1, nx-2, nx-2);
%–– build Y-implicit matrix for rows 2…ny−1
ey = ones(ny-2,1);
Ay = spdiags([ -r_y*ey, (1+2*r_y)*ey, -r_y*ey ], -1:1, ny-2, ny-2);

% while any(burning(:)) && n < nmax
while n < nmax
    
    % Increment time
    n = n + 1;
    temp_n = temp;

    Tmid = temp_n;
    for j = 2:ny-1
      RHS = temp_n(j,2:end-1) ...
          + r_y*( temp_n(j+1,2:end-1) - 2*temp_n(j,2:end-1) + temp_n(j-1,2:end-1) );
      Tmid(j,2:end-1) = Ax \ RHS';
    end
    
    % full-step: implicit in Y, explicit in X
    for i = 2:nx-1
      RHS = Tmid(2:end-1,i) ...
          + r_x*( Tmid(2:end-1,i+1) - 2*Tmid(2:end-1,i) + Tmid(2:end-1,i-1) );
      temp(2:end-1,i) = Ay \ RHS;
    end

    test_time_history(:,n) = [temp(time_history_j,time_history_i); burning(time_history_j,time_history_i); material(time_history_j,time_history_i)];
    
    % 1) cache the old burn‐map
    prev_status = node_status;

    % 2) find all cells hot enough to ignite
    grid_ignite = (temp >= temp_ign) & (material > 0);

    % 3) (optional) only allow ignition if touching the existing flame
    neigh = false(ny,nx);
    neigh(:,2:end)   = neigh(:,2:end)   | prev_status(:,1:end-1); % west
    neigh(:,1:end-1) = neigh(:,1:end-1) | prev_status(:,2:end);   % east
    neigh(2:end,:)   = neigh(2:end,:)   | prev_status(1:end-1,:); % south
    neigh(1:end-1,:) = neigh(1:end-1,:) | prev_status(2:end,:);   % north

    % 4) bona fide “just ignited” cells
    newIgn = grid_ignite & neigh & ~prev_status;

    % 5) update the burn map
    node_status = prev_status | newIgn;

    % 6) assign local WSB rates only to those newly‐ignited cells
    [Js,Is] = find(newIgn);
    for k=1:numel(Js)
      j = Js(k); i = Is(k);
      [rb_map(j,i),~,Tf_map(j,i)] = project_2_wsb_function(P, temp_n(j,i));
    end

    % 7) keep those burning at their flame-temperature
    temp(newIgn)  = Tf_map(newIgn);
    burning       = node_status & (material > 0);
    temp(burning) = Tf_map(burning);

    %node_status(material == 0) = 0;

    plot_refresh = 50;
    if uint16(n/plot_refresh) == n/plot_refresh % refresh the plot every 50 time steps to save time     
        clf
        
        subplot(1,3,1)
        contourf(X, Y, temp,  10, 'LineColor','none')
        title(sprintf('T @ %g s', n*dt))
        xlabel('x (cm)')
        ylabel('y (cm)')
        colorbar
        axis equal tight
        
        subplot(1,3,2)
        contourf(X, Y, material, 10, 'LineColor','none')
        title(sprintf('Thickness @ %g s', n*dt))
        xlabel('x (cm)')
        ylabel('y (cm)')
        axis equal tight
        
        subplot(1,3,3)
        contourf(X, Y, node_status, 2, 'LineColor','none')
        title(sprintf('Burning @ %g s', n*dt))
        xlabel('x (cm)')
        ylabel('y (cm)')
        axis equal tight
        
        pause(0.1)
    end
    

    % check for convergence
    stillSolid = ~node_status;  
    if ~any(burning(:))
        err = max(abs(temp(stillSolid) - temp_n(stillSolid)));
        if err <= tolerance
            fprintf("Converged at t=%g s after all burning ended (err=%g)\n", n*dt, err);
            break
        end
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

