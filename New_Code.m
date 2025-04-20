
% Numerically solve for transient conducton problem
clear all; close all; clc;

% obtain the input parameters
L = 0.1;                %Length of material [x-direction] (m)
W = 0.05;               %Width of material [y-direction] (m)
D = 0.05;               %Depth of material (m)
dx = 0.0025;            %Size of each element in x direction (m)
dy = 0.0025;            %Size of each element in y direction (m)
tmax = 300;              %Maximum iteration time (s)
dt = 0.01;              %Size of each time step (s)
epsilon = 10^(-4);      %Target Error to Exit Code Early
alp = 7.937*10^(-8);    %Thermal diffusivity of Condensed Phase (m^2/s)
rho = 1800;             %Density of Condensed Phase (kg/m^3)
E_c = 400*1000;         %Heat release of Condensed Phase (J/kg)
T_int = 25;             %Initial Starting Temperature (C)

alp_a = 7.937*10^(-8);    %Thermal diffusivity of air at 2000 C (m^2/s)


r_x = alp*dt/dx^2;
r_y = alp*dt/dy^2;
r_x_a = alp_a*dt/dx^2;
r_y_a = alp_a*dt/dy^2;
V = dx*dy*D;            %Calculates the volume of each element (m^3)         
m = rho*V;              %Calculates the mass of each element (kg)
E_t = E_c*m;            %Calculates total heat release of each element (J)

% create the x, y meshgrid based on dx, dy
nx = uint32(L/dx + 1);
ny = uint32(W/dy + 1);
[X,Y] = meshgrid(linspace(0,L,nx),linspace(0,W,ny));

% set initial and boundary conditions                
T = T_int*ones(ny,nx);    %Creates a matrix of each point at initial temperature
T(ny-1, 2) = 2000;        %Initial Starting points for burning surface temperatures
T(ny-1, 3) = 2000;
T(ny-1, 4) = 2000;
T(ny-2, 2) = 2000;
T(ny-2, 3) = 2000;
T(ny-2, 4) = 2000;
T(ny-3, 2) = 2000;
T(ny-3, 3) = 2000;
T(ny-3, 4) = 2000;

%Calculate the energy consumption rate of each element when burning
% r_b = 3.5*10^(-4);        %Steady-State Burn Rate (m^2/s)
r_b = 0.035;        %Steady-State Burn Rate (m^2/s)
T_burn = 2000;            %Steady-State burn temperature (C)
T_ig = 300;               %Ignition Temperature (C)
V_b = r_b*D;              %Volume Burn Rate (m^3/s)
t_b = V/V_b;              %Time for an element to fully burn (s)
E = t_b*ones(ny,nx);      %Creates a matrix of the stored energy of each point which may be consumed        


% iteration, march in time
n = 0; 
nmax = uint32(tmax/dt);
while n < nmax
    n = n + 1;

    for j = 1:ny
        for i = 1:nx
            if T(j,i)>= T_ig && E(j,i)>0        %Only burns if temperature is above ignition and there is energy
               T(j,i) = T_burn;
                E(j,i) = E(j,i)-dt;             %Subtract some of the energy which was used up
            
            end   
        end
        
    end

    
    
    T_n = T;

 
    for j = 2:ny-1
        for i = 2:nx-1
          
                if E(j-1,i)<= 0 && E(j+1,i)<= 0 && E(j,i-1)<= 0 && E(j,i+1)<= 0 && E(j,i)<=0
                T(j,i) = T_n(j,i) + r_x_a*(T_n(j,i+1)-2*T_n(j,i)+T_n(j,i-1))...
                + r_y_a*(T_n(j+1,i)-2*T_n(j,i)+T_n(j-1,i));
     
                else if E(j-1,i)<= 0 && E(j+1,i)<= 0 && E(j,i-1)<= 0 && E(j,i)<=0
                T(j,i) = T_n(j,i) + (r_x_a*T_n(j,i+1)-2*r_x*T_n(j,i)+r_x*T_n(j,i-1))...
                + r_y*(T_n(j+1,i)-2*T_n(j,i)+T_n(j-1,i));
           

                else if E(j-1,i)<= 0 && E(j+1,i)<= 0 && E(j,i+1)<= 0 && E(j,i)<=0
                T(j,i) = T_n(j,i) + (r_x_a*T_n(j,i+1)-r_x_a*2*T_n(j,i)+r_x*T_n(j,i-1))...
                + r_y_a*(T_n(j+1,i)-2*T_n(j,i)+T_n(j-1,i));

                else if E(j-1,i)<= 0 && E(j,i-1)<= 0 && E(j,i+1)<= 0 && E(j,i)<=0
                T(j,i) = T_n(j,i) + r_x_a*(T_n(j,i+1)-2*T_n(j,i)+T_n(j,i-1))...
                + (r_y*T_n(j+1,i)-2*r_y_a*T_n(j,i)+r_y_a*T_n(j-1,i));

                else if E(j+1,i)<= 0 && E(j,i-1)<= 0 && E(j,i+1)<= 0 && E(j,i)<=0
                T(j,i) = T_n(j,i) + r_x_a*(T_n(j,i+1)-2*T_n(j,i)+T_n(j,i-1))...
                + (r_y_a*T_n(j+1,i)-2*r_y_a*T_n(j,i)+r_y*T_n(j-1,i));

                else
                T(j,i) = T_n(j,i) + r_x*(T_n(j,i+1)-2*T_n(j,i)+T_n(j,i-1))...
                + r_y*(T_n(j+1,i)-2*T_n(j,i)+T_n(j-1,i));  

           

            end
            end
            end
            end
            
          end
           
         
        end

        
    end
    
    for i = 1:nx
        j = ny;
        T(j,i) = T(j-1, i);
    end

    for i = 1:nx
        j = 1;
        T(j,i) = T(j+1, i);
    end
   

    for j = 1:ny
        i = nx;
        T(j,i) = T(j, i-1);
    end
   
    for j = 1:ny
        i = 1;
        T(j,i) = T(j, i+1);
    end
   

    if uint16(n/50) == n/50 % refresh the plot every 50 time steps to save time     
        contourf(X,Y,T,10);
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


