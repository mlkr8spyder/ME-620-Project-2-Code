%% ME 620 Project 2: 2D Transient Heat Conduction Model for the Flame Spread of HMX
% Primary Author: Matthew Kelley
% Supporting Authors: Ethan Hund, Kaiden Nelson, Stephanie Farney
% Date: 4/26/25

clear all; close all; clc;

%% Solver Input parameters
grid_length = 4; % cm  
grid_height = 2; % cm  
dx = 0.05; % cm  (0.5 mm spacing)
dy = 0.05; % cm  (0.5 mm spacing)
time_max = 20;
dt = 0.0005;
n = 0; 
nmax = uint32(time_max/dt);
tolerance = 10^(-4);

%% Solid HMX Parameters
k = 0.2 * 1;
cp = 1400;
rho = 1800;
alpha = (k/(cp*rho)) * 1e4;
r_x = alpha*dt/dx^2;
r_y = alpha*dt/dy^2;

%% Gas HMX Parameters
k_gas = 0.07; % W/(m·K)
cp_gas = 1005;  % J/(kg·K)
rho_gas = 1.2; % kg/m^3
alpha_g = (k_gas/(cp_gas*rho_gas)) * 1e4;
r_x_g = alpha_g*dt/dx^2;
r_y_g = alpha_g*dt/dy^2;

%% WSB Initialize
P = 1; % Pressure [atm]
temp_ign = 500; % Ignition temp
[rb_wsb, Tf_test, temp_flame_wsb]   = project_2_wsb_function(P,temp_ign);
temp_flame = temp_flame_wsb;

%% Create the x, y meshgrid based on dx, dy
nx = uint32(grid_length/dx + 1);
ny = uint32(grid_height/dy + 1);

% Set the left and right buffer nodes to 10% of Y length
boundary_buffer = 0.1 * ny;
nx = nx + 2*boundary_buffer;
[X,Y] = meshgrid(linspace(0,grid_length,nx),linspace(0,grid_height,ny));

%% Material Map
thickness = 0.2; % thickness of material [cm]
material = ones(ny,nx).*thickness;

%% Initial and Boundary conditions
temp_int = 293;
temp = temp_int*ones(ny,nx);

% Logical Matrix to handle if node is active
node_status = zeros(ny,nx);

sides_inhibitor = 1; % Inhibitor BCs

% Initial Flame location from 48% to 52% along the Y length
boundary_flame_y_start = round(ny*0.48);
boundary_flame_y_end = round(ny*0.52);

material(1,:) = 0;    % top row
material(end,:) = 0;    % bottom row
material(:,1:boundary_buffer) = 0;    % left  column
material(:,end-(boundary_buffer-1):end) = 0;    % right column

% Set initial boundary to the ignition temperature
temp((boundary_flame_y_start):(boundary_flame_y_end),boundary_buffer+1) = temp_ign;

% Initialize the solid gas ratio, used for calculating the thermal
% diffusivity when there are two phases
solid_gas_ratio = material ./ thickness;

%% Burn Rate and Flame-Temp Maps
rb_map = zeros(ny,nx);
Tf_map = temp_int * ones(ny,nx);
t_surface_map = temp_int * ones(ny,nx);

%% Initial Ignition and Burning Logic Management
% Determines if a grid point should ignite by comparing temperature to the
% ingition temperature and making sure the node has material available
grid_ignite = (temp >= temp_ign) & (material > 0);

% Sets new node status based off previous node status and new grid ignition
node_status = node_status | grid_ignite;

% Is redudndant, but an artifact of a previous version of this code so it
% gets to stay in since burning is used in other areas
burning = node_status & (material > 0);

% Node status before iterating
node_status_n0 = node_status;

%% Intital burning calculations using WSB
% Intializes the BC's in the burning rate, surface temperature, and flame
% temperature maps
[j_list0, i_list0] = find(node_status_n0);
for k = 1:numel(i_list0)
    i_k0 = i_list0(k);
    j_k0 = j_list0(k);
    [rb_k0, t_surf_k0, t_flame_k0] = project_2_wsb_function(P,temp(j_k0,i_k0));
    
    rb_map(j_k0, i_k0) = rb_k0;
    Tf_map(j_k0, i_k0) = t_flame_k0;
    t_surface_map(j_k0,i_k0) = t_surf_k0;
end

%% Time Loop, runs until the max time is hit or there are no more burning nodes
while any(burning(:)) && n < nmax    
    %% Increment time
    n = n + 1;
    current_time = n*dt;
    temp_n = temp;

    %% Finite Difference Solver
    for j = 2:ny-1
        for i = 2:nx-1
            % If the node is burning, don't recalculate the nodes
            % temperature as it is fixed in this model
            if node_status(j,i) == 1
                temp(j,i) = temp_flame;
            else
                % Averaged thermal diffusivity for each node depending on
                % how much of the neighboring node is solid or gas
                rx_east = solid_gas_ratio(j,i+1)*r_x + (1-solid_gas_ratio(j,i+1))*r_x_g;
                rx_west = solid_gas_ratio(j,i-1)*r_x + (1-solid_gas_ratio(j,i-1))*r_x_g;
                ry_south = solid_gas_ratio(j+1,i)*r_x + (1-solid_gas_ratio(j+1,i))*r_x_g;
                ry_north = solid_gas_ratio(j-1,i)*r_x + (1-solid_gas_ratio(j-1,i))*r_x_g;
                
                % Node of interest and neighboring nodes
                temp_P = temp_n(j,i);
                temp_east = temp_n(j,i+1);
                temp_west = temp_n(j,i-1);
                temp_south = temp_n(j+1,i);
                temp_north = temp_n(j-1,i);
                
                % Used for debugging, not used in solver
                east_flux = rx_east*(temp_east-temp_P);
                west_flux = -rx_west*(temp_P-temp_west);
                south_flux = ry_south*(temp_south-temp_P);
                north_flux = -ry_north*(temp_P-temp_north);
                
                % Finite difference. 0 represents where we might add a heat
                % source term in later versions
                temp(j,i) = temp_P ...
                    + rx_east*(temp_east-temp_P) - rx_west*(temp_P-temp_west) ...
                    + ry_south*(temp_south-temp_P)   - ry_north*(temp_P-temp_north) ...
                    + 0;
            end
        end
    end

    %% Burning Node Management
    % Determines which nodes are currently burning
    burning = node_status & (material > 0);

    % Subtracts the material that burned away from during the time step
    % based upon that calculated burning rate for that node using WSB
    material = material - burning.*(rb_map.*dt);

    % Sets the minium material to 0 if the material burned away would make
    % it negative
    material(material < 0) = 0;

    % Calculate the ratio of solid to gas in the node
    solid_gas_ratio = material ./ thickness;

    % Turns off the nodes without any material left
    node_status (material == 0) = 0;
    
    % Save current node status to compare to new node status
    node_status_n0 = node_status;

    % Determine which new nodes ignite
    grid_ignite = (temp >= temp_ign) & (material > 0);

    % Creates new node status based off previous node status and new
    % ignition nodes
    node_status = node_status | grid_ignite;

    % Determines if the ignition of a node was new. Used to limit how many
    % times WSB is called, saving some computational time
    new_ignition = (node_status_n0 == 0) & (node_status == 1);

    %% Estimate the flame spread speed
    % Calculates the flame spread in cm/s speed by looking at how long it takes the
    % flame to reach 1 cm, 2 cm, 3 cm, and 4 cm
    if new_ignition(30, boundary_buffer+(((nx - 2*boundary_buffer)-1)/4)) == 1
        flame_speed_1 = 1 / (n*dt);
    elseif new_ignition(30, boundary_buffer+(2*((nx - 2*boundary_buffer)-1)/4)) == 1
        flame_speed_2 = 2 / (n*dt);
    elseif new_ignition(30, boundary_buffer+(3*((nx - 2*boundary_buffer)-1)/4)) == 1
        flame_speed_3 = 3 / (n*dt);
    elseif new_ignition(30, boundary_buffer+(4*((nx - 2*boundary_buffer)-1)/4)) == 1
        flame_speed_4 = 4 / (n*dt);
    end

    %% Calculate WSB for each new ignition
    [j_list, i_list] = find(new_ignition);
    % Runs through the new ignition nodes and saves the values from WSB to
    % a map. By establishing a map like this, we save a lot of computation
    % time.
    for k = 1:numel(i_list)
        i_k = i_list(k);
        j_k = j_list(k);
        [rb_k, t_surf_k, t_flame_k] = project_2_wsb_function(P,temp(j_k,i_k));

        rb_map(j_k, i_k) = rb_k;
        Tf_map(j_k, i_k) = t_flame_k;
        t_surface_map(j_k,i_k) = t_surf_k;
    end

    %% Plotting at defined time step intervals
    plot_refresh = 500;
    if uint16(n/plot_refresh) == n/plot_refresh % refresh the plot every plot_refresh time steps to save time     
        clf
        
        % We can plot the entire grid or just the propellant

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
end

%% Calculate the Average Flame Speed
average_flame_speed = mean([flame_speed_1 flame_speed_2 flame_speed_3 flame_speed_4]);
fprintf('The average flame speed is %f cm/sec',average_flame_speed)

