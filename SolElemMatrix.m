function [elematSol] = SolElemMatrix(f,eID,mesh)
% Solution Matrix
%This function generates a local element vector for the solution vector
% for an arbitrary element Ne between points xmin and xmax

J = abs(mesh.elem(eID).J);        % Returns Jacobian value from mesh data structure

% Initialises local element matrix
Int = zeros(3,1);

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
        
        Int(n) = gq.gsw(i) * (f * psi(n) * J) + Int(n);
        
    end
    
end

elematSol = Int;

end