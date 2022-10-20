function [psi] = EvalBasis(lnid,xipt)
% EvalBasis returns values of basis functions
% Returns the values of the quadratic basis functions for a
% specificed local node id (lnid = 0 or 1 or 2) and xipt

switch lnid
    case 1
psi = xipt.*(xipt-1)/2;
    case 2
psi = 1 - xipt.^2;
    case 3
psi = xipt.*(xipt+1)/2;
end
end


