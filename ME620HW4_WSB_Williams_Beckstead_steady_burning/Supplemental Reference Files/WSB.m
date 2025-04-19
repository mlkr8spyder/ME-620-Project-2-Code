% Original Author: Andrew Rettenmaier
% Oct. 10th, 2011
% ME 697C

% Version 2 Author: Allison Murray
% 19 Oct 2016

% Problem 3, Whitson parameters

% Compute deflagration using WSB model and compares to WDB (williams 1973)

%% User Inputs
P = 1;                    % atm
T0 = 293:20:673;           % K, desired initial temperature
%T0 = 298;


%% Constants
for kk = 1:length(T0)
    for ii = 1:length(P)
        %% Model Inputs
        % for ZN stability fit from Whitson et al
        D.Qg    = 3827*1000;                                     % J/kg,        gas phase heat release
        D.Qc    = -109*1000;                                      % J/kg,        condensed phase heat release
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
        
        %% Dennison & Baum        
        md_DB(ii,kk) = sqrt(2*D.kg*D.BgW*D.MW^2*D.p^2*D.cp*D.Tf^4/(D.Eg^2*D.Qg^2))...% kg/m^2-s, mass burning rate
            *exp(-D.Eg/(2*D.R*D.Tf));
        rd_DB(ii,kk) = md_DB(ii,kk)/D.rhoc;                                         % m/s, burning rate
        rd_DB(ii,kk) = rd_DB(ii,kk)*100;                                     % cm/s, burning rate
        Ts_DB(ii,kk) = (-D.Ec/D.R)/(log(md_DB(ii,kk)/(D.Ac*D.rhoc)));
        xgd_DB(ii,kk) = D.kg*(md_DB(ii,kk)*D.cp*log(D.Qg/(D.cp*(Ts_DB(ii,kk)-D.T0)-D.Qc)))^-1;
    end
end

%% Plotting
figure(1);
T0_cels = T0-273;
plot(T0_cels,rd,'r',T0_cels,rd_DB,'b');
title('ME 697C, HW 3, Problem 3, RDX - Allison Murray');
ylabel('Burning rate, cm/s');
xlabel('Initial Temperature, C');
legend('WSB','WDB','location','best');
axis tight
figure(2);
plot(T0_cels,xgd,'r',T0_cels,xgd_DB,'b');
title('ME 697C, HW 3, Problem 3, RDX - Allison Murray');
ylabel('Flame Standoff, mm');
xlabel('Initial Temperature, C');
legend('WSB','WDB','location','best');
axis tight

% 
% 
% 
