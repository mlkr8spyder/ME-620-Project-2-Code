 
% Andrew Rettenmaier
% Oct. 10th, 2011
% ME 697C
% P  = [1 5 10 20 70 80];
% T0 = 290:2:310;
select = 1; %select = 1 for HMX select = 2 %for AP
%% Inputs
if select == 1
    % HMX
    A.As = 0.9e10;        % g/cm^2-s
    A.Es = 50000;         % cal/mol
    A.Qc = -200;           % cal/g
    A.cp = 0.3;           % cal/g-K
    A.rho = 1.91;         % g/cm^3
    A.lamb = 0.0003;      % cal/cm-K-s
    A.delt = 1.6;         % ??
    A.R    = 1.985;       % cal/mol-K
    A.T0   = 298;
    r0   = 1.09;          % cm/s
    A.Tf   = 3350;        % K
    P0      = 68*101325*10;         % kg-m/s^2-m^2
    A.Qf   = A.cp*(A.Tf - A.T0) - A.Qc;
else
    % AP
    A.As = 0.3e6;         % g/cm^3-s
    A.Es = 22000;         % cal/mol
    A.Qc = -120;          % cal/g
    A.cp = 0.3;           % cal/g-K
    A.rho  = 1.95;        % g/cm^3
    A.lamb = 0.0003;      % cal/cm-K-s
    A.delt = 1.8;         % ??
    A.R    = 1.985;       % cal/mol-K
    r0     = 0.83;        % cm/s
    A.Tf   = 1350;        % K
    A.T0   = 298;         % K
    P0     = 68*101325*10;  % Pa
    A.Qf   = A.cp*(A.Tf - A.T0) - A.Qc;
end

% solve at reference condition
Ts0f = @(y)(log(r0*A.rho*sqrt(y - A.T0)/A.As) + A.Es/(2*A.R*y));
options = optimset('MaxFunEvals',10000,'TolX',1e-5,'TolFun',1e-5,'MaxIter',10000);
y0 = -A.Es/(A.R*log(r0*A.rho/A.As));

%Ts0   = fsolve(Ts0f,y0,options);
Ts0   = y0;
chi = -log((A.cp*(Ts0-A.T0) - A.Qc)/A.Qf);
A.k   = A.cp*(r0*A.rho)^2/(chi*A.lamb*P0^A.delt);  % Assume this is now constant

%mtest = (A.As/(sqrt(Ts0 - A.T0)))*exp(-A.Es/(2*A.R*Ts0)); % g/cm^2-s should match reference mass burning rate
mtest = (A.As*exp(-A.Es/(A.R*Ts0))); % g/cm^2-s should match reference mass burning rate

for jj = 1:length(T0)
    for ii = 1:length(P)
        A.T0   = T0(jj);         % K
        A.P    = P(ii)*101325*10;   % Pa
        
        %% Solve at desired pressure
        x0 = [.6*A.rho 700 chi];
        becksteadf = @(x)becksolve(x,A);
        options = optimset('MaxFunEvals',10000,'TolX',1e-5,'TolFun',1e-5,'MaxIter',10000);
        
        [x,fval] = fsolve(becksteadf,x0,options);
        
        %% Output
        mb(ii,jj) = (x(1)); % g/cm^2-s
        rb(ii,jj) = (x(1))/A.rho;       % cm/s
        Tsb(ii,jj) = (x(2));      % K
        chib(ii,jj) = (x(3));
        xf(ii,jj)    = mb(ii,jj)/(A.k*A.P^A.delt);
    end
end



