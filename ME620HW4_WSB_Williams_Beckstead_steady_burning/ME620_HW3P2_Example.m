% Original Author: Andrew Rettenmaier
% Oct. 10th, 2011
% ME 697C
% Version 2 Author: Allison Murray
% 19 Oct 2016
% Cleaned up by Jim Plotzke
% February 2024
clear, clc, clf, close all

% Problem 3, Whitson parameters
% Compute deflagration using WSB model and compares to WDB (williams 1973)

%% User Inputs
%##################################################################################
% HEY YOU, CHANGE THESE VALUES TO WHAT THEY NEED TO BE
P = 50; % atm
T0 = 300:1:400; %K
% P.S., you'll need the Optimization Toolbox
% If you don't know how to download it, click on the HOME tab, then look on the ribbon for "Add-Ons" -> "Get Add-Ons"

%####################################################################################
%T0 = 673;           % K, desired initial temperature
%% q_bar = 200;                 %W/cm^2
% D.rhoc  = 1800/(100^3); 
% r_b = 0.165;
% D.cp    = 1.4*1000;   
% 
% T0 = 298 + q_bar/(D.rhoc*r_b*D.cp);

%% Constants
for kk = 1:length(T0)
    for ii = 1:length(P)
        %% Model Inputs
        % for ZN stability fit from Whitson et al
        D.Qg    = 3018*1000;                                     % J/kg,        gas phase heat release
        D.Qc    = 400*1000;                                      % J/kg,        condensed phase heat release
        D.cp    = 1.4*1000;                                      % J/kg-K,      condensed phase specific heat
        D.kg    = 0.07;                                          % W/m-K,       gas phase conductivity
        D.Ac    = 1.637e15;                                      % 1/s,         Arrhenius pre-factor
        D.kc    = 0.2;                                           % W/m-K,       condensed phase conductivity
        D.Ec    = 176*1000;                                      % J/mol,       condensed phase activation energy
        D.Eg    = 167.36*1000;                                      % J/mol,       gas phasae activiation energy (for williams model)
        D.rhoc  = 1800;                                          % kg/m^3,      condensed phase density
        D.MW    = 34.2/1000;                                     % kg/mol,     condensed phase molecular weight
        D.Bg    = 1.6e-3;                                        % m^3/kg-s-K^2,Arrhenius pre-factor
        D.BgW   = 377;                                           % m^3/kg-s-K^2,Arrhenius pre-factor for williams model
        D.T0    = T0(kk);                                        % K,           Initial temperature (specify above)
        D.Tf    = D.T0+(D.Qc + D.Qg)/D.cp;                       % K,           Flame temperature
        D.R     = 8.314;                                         % J/mol-K,     Universal gas constant
        D.mr    = 1;                                             % m-kg/s-m^3,  reference mass burning rate
        D.p     = P(ii)*101325;                                  % Pa,          burning pressure
        %% Non-dimensional values
        A.Ec   = D.Ec/(D.R*(D.Tf-D.T0));                         % Condensed phase activation energy
        A.Eg   = D.Eg/(D.R*(D.Tf-D.T0));                         % Gas phase activation energy
        A.Qc   = D.Qc/(D.cp*(D.Tf-D.T0));                        % Condensed phase heat release
        A.Qg   = D.Qg/(D.cp*(D.Tf-D.T0));                        % Gas phase heat release
        A.Ac   = D.kc*D.rhoc*D.Ac/(D.cp*D.mr^2);                 % Pre-factor
        A.Dg   = D.kg*D.Bg*D.p^2*D.MW^2/(D.mr^2*D.R^2*D.cp);     % Damkoehler number
        A.Tf   = D.Tf/(D.Tf-D.T0);                               % Flame temperature
        A.T0   = D.T0/(D.Tf-D.T0);                               % Initial temperature
        
%% WSB
        Tsguess = 0.3491;                                        % Guess at nd surface temperature
        mguess  = sqrt((A.Ac*Tsguess^2*exp(-A.Ec/Tsguess)/...    % mass burning rate at Tsguess
            (A.Ec*(Tsguess-A.T0-A.Qc/2))));
        xgguess = 2/(sqrt(mguess^2 + 4*A.Dg) - mguess);          % flame stand-off distance at Tsguess
        
        x0 = [Tsguess mguess xgguess];                           % initial guess vector
        
        WSBf = @(x)WSBsolve(x,A);                                % creates function handle for WSBsolve
        options = optimset('MaxFunEvals',10000,'TolX',...        % options for fsolve function
            1e-10,'TolFun',1e-10,'MaxIter',10000);
        [x,fval] = fsolve(WSBf,x0,options);                      % solve non-linear coupled equations
        
        %% Create dimensional forms
        Ts = x(1);                               % Surface temperature
        m  = x(2);                               % mass burning rate
        xg = x(3);                               % flame stand-off distance
        
        md(ii,kk)  = m*D.mr*1000/100^2;           % g/cm^2-s,  mass burning rate (6.82 kg/cm^2-s)
        rd(ii,kk)  = m*D.mr*100/D.rhoc;           % cm/s,      burning rate
        Tsd(ii,kk) = Ts*(D.Tf-D.T0);              % K,         Surface temperature
        xgd(ii,kk) = xg*(D.kg/(D.mr*D.cp))*1000;  % mm,        stand off distance
        D.Ti = D.T0 + D.Qc/D.cp + D.Qg/(D.cp*(md(ii,kk)*xgd(ii,kk)*D.cp/D.kg+1));
        sigma_p(ii,kk) = 0.5*(D.rhoc*D.Ec*(D.cp*(D.Ti-D.T0)-D.Qc/2))^(5/4)*((D.rhoc*D.Ec*... %% Temperature Sensitivity
            (D.cp*(D.Ti-D.T0)-D.Qc/2))*(2*D.Ac*D.R*D.T0*D.kc*exp(-1*D.Ec/(D.R*D.Ti)))+...
            D.Ac*D.R*D.T0^2*D.kc*exp(-1*D.Ec/(D.R*D.Ti))*D.rhoc*D.Ec*D.cp);
        
        q_crit(ii,kk) = D.kc/100*(Tsd(ii,kk)-D.T0)/rd(ii,kk);
        
        
        %% Dennison & Baum        
        md_DB(ii,kk) = sqrt(2*D.kg*D.BgW*D.MW^2*D.p^2*D.cp*D.Tf^4/(D.Eg^2*D.Qg^2))...% kg/m^2-s, mass burning rate
            *exp(-D.Eg/(2*D.R*D.Tf));
        rd_DB(ii,kk) = md_DB(ii,kk)/D.rhoc;                                         % m/s, burning rate
        rd_DB(ii,kk) = rd_DB(ii,kk)*100;                                     % cm/s, burning rate
        Ts_DB(ii,kk) = (-D.Ec/D.R)/(log(md_DB(ii,kk)/(D.Ac*D.rhoc)));
        xgd_DB(ii,kk) = D.kg*(md_DB(ii,kk)*D.cp*log(D.Qg/(D.cp*(Ts_DB(ii,kk)-D.T0)-D.Qc)))^-1;
    end
end

%% Sensitivity Analysis
% sens_wsb = (diff(xgd)./diff(P)')./xgd(1:length(xgd)-1);
% sens_wdb = (diff(xgd_DB)./diff(P)')./xgd_DB(1:length(xgd_DB)-1);
% sens1_wsb = sens_wsb(1);
% senslast_wsb = sens_wsb(length(sens_wsb));
% sens1_wdb = sens_wdb(1);
% senslast_wdb = sens_wdb(length(sens_wdb));
% r = diff(Tsd)./diff(T0);
% k = (Tsd(1:length(Tsd)-1)-T0(1:length(T0)-1)).*(diff(log(md))./diff(T0));
% r_DB = diff(Ts_DB)./diff(T0);
% k_DB = (Ts_DB(1:length(Ts_DB)-1)-T0(1:length(T0)-1)).*(diff(log(md_DB))./diff(T0));
%% Plotting
%#############################################################################################################
% FOR PROBLEM 3, YOU'LL WANT A SCREENSHOT OF THIS PLOT AT THE TWO DIFFERENT STARTING PRESSURES
figure(1);
hold on
T0_cels = T0-273.15;
plot(T0_cels,rd,'r')
plot(T0_cels,rd_DB,'b')
plot(T0_cels, 0) %delete this
plot(T0_cels, 0) %delete this
hold off
ylabel('Burning rate, cm/s');
xlabel('Initial Temperature, C');
%edit the legend and title
title('Read through the code, find this plot, and edit it.');
legend('WSB','WDB','Why do you have to edit the legend and title?','So I know you ran the code instead of just copying the solution.','location','best');
axis tight
%#################################################################################################################

% IGNORE THE REST OF THESE PLOTS FOR PROBLEM 3
% figure(2);
% plot(T0_cels,xgd,'r',T0_cels,xgd_DB,'b');
% title('ME 697C, HW 3, Problem 3, HMX - Allison Murray');
% ylabel('Flame Standoff, mm');
% xlabel('Initial Temperature, C');
% legend('WSB','WDB','location','best');
% axis tight
% figure(1);
% plot(P, q_crit);
% title('ME 697C, HW 4, Problem 4, HMX - Allison Murray');
% ylabel('Critical Energy (J/cm^2)');
% xlabel('Pressure (atm)');
% figure(1);
% plot(T0(1:length(T0)-1),r);
% title('ME 697C, HW 4, Problem 6, HMX, p=0.5 atm - Allison Murray');
% ylabel('r - ZN sensitivity');
% xlabel('Temperature(K)');
% figure(2);
% plot(T0(1:length(T0)-1),k);
% title('ME 697C, HW 4, Problem 6, HMX, p=0.5 am - Allison Murray');
% ylabel('k - Temperature sensitivity');
% xlabel('Temperature(K)');
% figure(3);
% plot(k,r);
% title('Stability plot, HW 4, Problem 6, p=0.5 atm - Allison Murray');
% ylabel('ZN sensitivity');
% xlabel('Temperature sensitivity');
% figure(4);
% plot(T0(1:length(T0)-1),r_DB);
% title('ME 697C, HW 4, Problem 6, HMX (WDB), p=0.5 atm - Allison Murray');
% ylabel('r - ZN sensitivity');
% xlabel('Temperature(K)');
% figure(5);
% plot(T0(1:length(T0)-1),k_DB);
% title('ME 697C, HW 4, Problem 6, HMX (WDB), p=0.5 am - Allison Murray');
% ylabel('k - Temperature sensitivity');
% xlabel('Temperature(K)');
% figure(6);
% plot(k_DB,r_DB);
% title('Stability plot, HW 4, Problem 6, p=0.5 atm  (WDB)- Allison Murray');
% ylabel('ZN sensitivity');
% xlabel('Temperature sensitivity');
% 
% 
% 
