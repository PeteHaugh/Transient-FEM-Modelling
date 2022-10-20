function [elematLin] = LinElemMatrix(L,mesh)
% Linear Element Matrix
%This function generates a local 3-by-3 element matrix for the linear
%reaction operator for an arbitrary element Ne between points xmin and xmax 

J = mesh.elem.J;        % Returns Jacobian value from mesh data structure

% Initialises local element matrix
Int = zeros(3,3);

N=3;
[gq] = CreateGQScheme(N);

for i = 1:N
    
xipt = gq.xipts(i);

[psi] = [EvalBasis(1,xipt) EvalBasis(2,xipt) EvalBasis(3,xipt)];

for n = 1:3
    for m = 1:3
        
Int(n,m) = gq.gsw(i)* L * psi(n) * psi(m) * J + Int(n,m);

    end
end

end

elematLin = Int;

end

