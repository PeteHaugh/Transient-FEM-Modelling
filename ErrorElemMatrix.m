function [E] = ErrorElemMatrix(xvec, C, eID, t)
% L2 local error calculator
% This function calculates the L2 error between the numeric basis function solution 

% Initialise Gaussian quadrature scheme   
N=3;
[gq] = CreateGQScheme(N);

% Calculate psi for xi pt 1 and 3 (the middle point is negated)
psi1 = EvalBasis(1, gq.xipts);
psi2 = EvalBasis(3, gq.xipts);

% Determine local element values
x1 = xvec(eID);
x2 = xvec(eID + 1);

% Calculate anayltical solution for local element values
CE = TransientAnalyticSoln((x1 * psi1)+(x2 * psi2), t);

% Calculate numeric solution
C = (C(((2 * eID) - 1), t) * psi(1)) + (C(((2 * eID) + 1), t) * psi(2));

% Calculate error
E = (CE - C).^2;

end

