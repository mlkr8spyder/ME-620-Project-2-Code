% Numerically solve for transient conducton problem
clear all; close all; clc;

% Input parameters
grid_length = 4;    % cm  
grid_height = 2;     % cm  
dx = 0.05;            % cm  (0.5 mm spacing)
dy = 0.05;            % cm  (0.5 mm spacing)
time_max = 10;
dt = 0.0005;
n = 0; 
nmax = uint32(time_max/dt);
tolerance = 10^(-4);

%% Solid Parameters
k = 0.2 * 10;
cp = 1400;
rho = 1800;
alpha   = (k/(cp*rho)) * 1e4;
r_x     = alpha*dt/dx^2;
r_y     = alpha*dt/dy^2;

%% Gas Parameters
k_gas   = 0.07;       % W/(m·K)
cp_gas  = 1005*2;        % J/(kg·K)
rho_gas = 1.2;         % kg/m^3
alpha_g = (k_gas/(cp_gas*rho_gas)) * 1e4;
r_x_g   = alpha_g*dt/dx^2;
r_y_g   = alpha_g*dt/dy^2;

%% WSB Initialize
P = 1; % Pressure [atm]
temp_ign = 500; % Arbitrary ignition temp
[rb_wsb, Tf_test, temp_flame_wsb]   = project_2_wsb_function(P,temp_ign);
temp_flame = temp_flame_wsb;

%% Create the x, y meshgrid based on dx, dy
nx = uint32(grid_length/dx + 1);
ny = uint32(grid_height/dy + 1);
[X,Y] = meshgrid(linspace(0,grid_length,nx),linspace(0,grid_height,ny));

%% Material Map
thickness = 0.2; % thickness of material [cm]
material = ones(ny,nx).*thickness;

%% Initial and boundary conditions
temp_int = 293;
temp = temp_int*ones(ny,nx);

% Logical Matrix
node_status = zeros(ny,nx);

sides_inhibitor = 1;

if sides_inhibitor == 1
    boundary_buffer = 0.1 * ny;
    boundary_flame_y_start = 2;
    boundary_flame_y_end = ny-1;

    material(1,:) = 0;    % top row
    material(end,:) = 0;    % bottom row
    material(:,1:boundary_buffer) = 0;    % left  column
    material(:,end-(boundary_buffer-1):end) = 0;    % right column

    temp((boundary_flame_y_start):(boundary_flame_y_end),boundary_buffer+1) = temp_flame;
else 
    boundary_buffer = 0.1*ny;
    boundary_flame_y_start = boundary_buffer + 0.1 * (ny-(2*boundary_buffer));
    boundary_flame_y_end = boundary_buffer + 0.9 * (ny-(2*boundary_buffer));

    material(1:boundary_buffer,  :) = 0;    % bottom row
    material(end-(boundary_buffer-1):end,:) = 0;    % top row
    material(:,  1:boundary_buffer) = 0;    % left  column
    material(:,end-(boundary_buffer-1):end) = 0;    % right column

    temp((boundary_flame_y_start):(boundary_flame_y_end),boundary_buffer+1) = temp_flame;
end

%% Burn Rate and Flame-Temp Maps
rb_map = zeros(ny,nx);
Tf_map = temp_int * ones(ny,nx);
t_surface_map = temp_int * ones(ny,nx);

%% Initial Ignition
grid_ignite = (temp >= temp_ign) & (material > 0);
node_status = node_status | grid_ignite;
burning = node_status & (material > 0);
node_status_n0 = node_status;

test_time_history = zeros(3,nmax);
time_history_j = round(ny/2);
time_history_i = round(nx/2);

while any(burning(:)) && n < nmax
% while n < nmax
    
    % Increment time
    n = n + 1;
    current_time = n*dt;
    temp_n = temp;

    for j = 2:ny-1
        for i = 2:nx-1
            if burning(j,i) == 1
                temp(j,i) = temp_flame;
            else
                isGas = (temp_n(j,i) >= temp_ign) || (material(j,i) == 0);
                
                % pick diffusivities for each neighbor direction
                if isGas || material(j,i+1)==0, rx_east = r_x_g; else rx_east = r_x; end
                if isGas || material(j,i-1)==0, rx_west = r_x_g; else rx_west = r_x; end
                if isGas || material(j+1,i)==0, ry_south = r_y_g; else ry_south = r_y; end
                if isGas || material(j-1,i)==0, ry_north = r_y_g; else ry_north = r_y; end
                
                temp_P = temp_n(j,i);
                temp_east = temp_n(j,i+1);
                temp_west = temp_n(j,i-1);
                temp_south = temp_n(j+1,i);
                temp_north = temp_n(j-1,i);
    
                east_flux = rx_east*(temp_east-temp_P);
                west_flux = -rx_west*(temp_P-temp_west);
                south_flux = ry_south*(temp_south-temp_P);
                north_flux = -ry_north*(temp_P-temp_north);
    
                temp(j,i) = temp_P ...
                    + rx_east*(temp_east-temp_P) - rx_west*(temp_P-temp_west) ...
                    + ry_south*(temp_south-temp_P)   - ry_north*(temp_P-temp_north) ...
                    + 0;
            end
        end
    end

    test_time_history(:,n) = [temp(time_history_j,time_history_i); burning(time_history_j,time_history_i); material(time_history_j,time_history_i)];
    
    % Turn off nodes that have burnt out here
    burning = node_status & (material > 0);
    material = material - burning.*(rb_wsb*dt);
    material(material < 0) = 0;

    node_status (material == 0) = 0;
    
    % Search for values exceeding the ignition temperature
    node_status_n0 = node_status;
    grid_ignite = (temp >= temp_ign) & (material > 0);
    node_status = node_status | grid_ignite;
    new_ignition = (node_status_n0 == 0) & (node_status == 1);

    %Calculate WSB for each new ignition
    [j_list, i_list] = find(new_ignition);
    for k = 1:numel(i_list)
        i_k = i_list(k);
        j_k = j_list(k);
        [rb_k, t_surf_k, t_flame_k] = project_2_wsb_function(P,temp(j_k,i_k));

        rb_map(j_k, i_k) = rb_k;
        Tf_map(j_k, i_k) = t_flame_k;
        t_surface_map(j_k,i_k) = t_surf_k;
    end

    % Set any values exceeding ignition temperature to flame temp
    temp(node_status) = temp_flame;

    plot_refresh = 50;
    if uint16(n/plot_refresh) == n/plot_refresh % refresh the plot every 50 time steps to save time     
        clf
        
        % subplot(1,3,1)
        % contourf(X, Y, temp,  10, 'LineColor','none')
        % title(sprintf('T @ %g s', n*dt))
        % xlabel('x (cm)')
        % ylabel('y (cm)')
        % colorbar
        % axis equal tight
        propellant_only_x = [boundary_buffer + 1:nx-boundary_buffer];
        propellant_only_y = [boundary_buffer + 1:ny-boundary_buffer];
        subplot(1,3,1)
        contourf(X(propellant_only_y,propellant_only_x), Y(propellant_only_y,propellant_only_x), temp(propellant_only_y,propellant_only_x),  10, 'LineColor','none')
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
    % stillSolid = ~node_status;  
    % if any(stillSolid(:))
    %     err = max(abs(temp(stillSolid) - temp_n(stillSolid)));
    %     if err <= tolerance
    %          fprintf("Stopping at step %d (t=%g s) due to convergence: err=%g ≤ tol\n", ...
    %                 n, n*dt, err);
    %         break
    %     end
    % end
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

