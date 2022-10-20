function [ dpsidxi ] = EvalBasisGrad(lnid,xipt)
% EvalBasisGrad Returns gradient of basis functions
% Returns the gradients of the quadratic basis functions for a
% specificed local node id (lnid = 0 or 1 or 2) and xipt 
    
switch lnid
    case 1
dpsidxi = xipt-0.5;
    case 2
dpsidxi = - 2 * xipt;
    case 3
dpsidxi = xipt+0.5;
end

end

