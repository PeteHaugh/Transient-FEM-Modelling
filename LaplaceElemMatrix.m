function [elematLap] = LaplaceElemMatrix(D,eID,mesh)
% Laplace Element Matrix
%This function generates a local 3-by-3 element matrix for the diffusion
%operator for an arbitrary element Ne between points xmin and xmax

J = abs(mesh.elem(eID).J);        % Returns Jacobian value from mesh data structure
x = mesh.elem(eID).x;             % Returns local x values from mesh data structure
dXi = abs(2/(x(1) - x(2)));

% Initialises local element matrix
Int = zeros(3,3);

% Initialises Gaussian quadrature scheme
N=3;
[gq] = CreateGQScheme(N);

% Runs Gaussian quadrature scheme
for i = 1:N
    
    % Returns Xi points
    xipt = gq.xipts(i);
    
    % Returns the gradients of basis function
    [dpsidxi] = [EvalBasisGrad(1,xipt) EvalBasisGrad(2,xipt) EvalBasisGrad(3,xipt)];
    
    % Sums Gaussian quadrature values
    for n = 1:3
        for m = 1:3
            
            Int(n,m) = gq.gsw(i) * (D * dpsidxi(n) * dXi * dpsidxi(m) * dXi * J) + Int(n,m);
            
        end
    end
    
end

elematLap = Int;


end

