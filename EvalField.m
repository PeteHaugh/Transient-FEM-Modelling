function [nval] = EvalField(msh,field,eID,xipt)


 psi = [EvalBasis(0,xipt) EvalBasis(1,xipt)]; %Get vector containing basis function values at xipt
 vcoeff = [field(msh.elem(eID).n(1)); field(msh.elem(eID).n(2))]; %Get nodal values for element eID
 nval = psi*vcoeff; %Dot product the two together to sum up the coefficients*basis functions to give the final interpolated value

end


