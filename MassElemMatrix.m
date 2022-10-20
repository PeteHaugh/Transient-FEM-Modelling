function [M] = MassElemMatrix(eID,mesh)
% Mass Element Matrix
% This function generates a local 3-by-3 element matrix for the mass matrix
% for an arbitrary element Ne between points xmin and xmax 

J = mesh.elem(eID).J;        % Returns Jacobian value from mesh data structure

% Initialises local element matrix
Int = zeros(3,3);

% Initialises Gaussian quadrature scheme
N=3;
[gq] = CreateGQScheme(N);

% Runs Gaussian quadrature scheme
for i = 1:N
    
  % Returns Xi points   
xipt = gq.xipts(i);

% Returns the values of basis function
[psi] = [EvalBasis(1,xipt) EvalBasis(2,xipt) EvalBasis(3,xipt)];

% Sums Gaussian quadrature values
for n = 1:3
    for m = 1:3
        
Int(n,m) = gq.gsw(i) * psi(n) * psi(m) * J + Int(n,m);

    end
end

end

M = Int;


end

