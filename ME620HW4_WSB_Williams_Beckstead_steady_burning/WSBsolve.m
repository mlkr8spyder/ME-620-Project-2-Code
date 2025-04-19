function F = WSBsolve(x,A)

Ts   = x(1);
m    = x(2);
xg   = x(3);

F = [m^2*(A.Ec*(Ts - A.T0 - A.Qc/2)) - (A.Ac*Ts^2*exp(-A.Ec/Ts));
    (Ts - A.T0 - A.Qc)*(xg*m + 1) - A.Qg;
    xg*(sqrt(m^2 + 4*A.Dg) - m) - 2];
