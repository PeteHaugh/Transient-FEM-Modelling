function [elematLin] = LinOpElemMatrix(L,eID,mesh)
%% Linear Operator Element Matrix
%This function generates a local 2-by-2 element matrix for the linear
%reaction perator for an arbitrary element Ne between points xmin and xmax 

x = mesh.elem(eID).x;   % Returns local x values from mesh data structure

J = abs((x(2) - x(1))/2); % Calculates Jacobian from x values

%Int(n,m)


Int00 = L*J*2/3;
Int01 = L*J*1/3;
Int10 = Int01;
Int11 = L*J*2/3;

% Stores calculated values in element matrix
elematLin = zeros(2,2);
elematLin(1,1) = Int00;
elematLin(1,2) = Int01;
elematLin(2,1) = Int10;
elematLin(2,2) = Int11;

end

