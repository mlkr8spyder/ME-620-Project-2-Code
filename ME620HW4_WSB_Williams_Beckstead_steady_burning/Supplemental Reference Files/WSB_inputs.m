% Andrew Rettenmaier
% WSB INPUTS
% October 11th, 2011

% D.T0 and D.p may be vectors
% Current parameters for HMX

% for Pressure Response data (best matched to data from paper)
D.Qg    = 3171.5*1000;                                     % J/kg,        gas phase heat release
D.Qc    = 251.04*1000;                                      % J/kg,        condensed phase heat release
% for ZN stability fit from Ward et al
%D.Qg    = 3018*1000;                                     % J/kg,        gas phase heat release
%D.Qc    = 400*1000;                                      % J/kg,        condensed phase heat release
D.cp    = 1.4*1000;                                      % J/kg-K,      condensed phase specific heat
D.kg    = 0.07;                                          % W/m-K,       gas phase conductivity
D.Ac    = 1.637e15;                                      % 1/s,         Arrhenius pre-factor
D.kc    = 0.2;                                           % W/m-K,       condensed phase conductivity
D.Ec    = 176*1000;                                      % J/mol,       condensed phase activation energy
D.Eg    = 176*1000;                                      % J/mol,       gas phasae activiation energy (for williams model)
D.rhoc  = 1800;                                          % kg/m^3,      condensed phase density
D.MW    = 34.2/1000;                                     % kg/mol,     condensed phase molecular weight
D.Bg    = 1.6e-3;                                        % m^3/kg-s-K^2,Arrhenius pre-factor
D.BgW   = 436;                                           % m^3/kg-s-K^2,Arrhenius pre-factor for williams model
D.T0    = T0(kk);                                        % K,           Initial temperature (specify above)
D.Tf    = D.T0+(D.Qc + D.Qg)/D.cp;                       % K,           Flame temperature
D.R     = 8.314;                                         % J/mol-K,     Universal gas constant
D.mr    = 1;                                             % m-kg/s-m^3,  reference mass burning rate
D.p     = P(ii)*101325;                                  % Pa,          burning pressure