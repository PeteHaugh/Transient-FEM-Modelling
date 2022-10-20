function [GM,GV] = Dirichlet(No,GM,GV,BC)
% Dirichlet Boundary Conditions
% This function applies dirichlet boundary conditions to the values of the global matrix and global vector

% Applies boundary conditions to global vector
GV(1) = BC(1);
GV(end) = BC(2);

% Applies boundary conditions to global matrix
GM(1,:) = 0;
GM(No,:) = 0;
GM(1,1) = 1;
GM(No,No) = 1;

end
