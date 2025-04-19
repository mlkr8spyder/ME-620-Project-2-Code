function F = becksolve(x,A)
% Andrew Rettenmaier
% Oct. 10th, 2011
% ME 697C

m = x(1);
Ts = x(2);
chi = x(3);
% 
% F = [(log(m*sqrt(Ts - A.T0)/A.As) + A.Es/(2*A.R*Ts));
%      log((A.cp*(Ts - A.T0) - A.Qc)/A.Qf) + chi;
%      chi - A.cp*(m)^2/(A.lamb*A.k*A.P^A.delt)];

 F = [m - A.As*exp(-A.Es/(A.R*Ts));
     (A.cp*(Ts - A.T0) - A.Qc) - A.Qf*exp(-chi);
     chi - A.cp*(m)^2/(A.lamb*A.k*A.P^A.delt)];


 

