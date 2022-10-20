function [GM,mesh,S] = GMGen(Ne,D,L,xmin,xmax,f,mode)
%% Global Matrix Generator
% This function creates a global matrix (GM) and returns the global
% matrix, the mesh data structure and the solution vector (S)

mesh = OneDimLinearMeshGen(xmin,xmax,Ne);
No = Ne + 1;        % Converts No. of elements to No. of nodes
GM = zeros(No);     % CInitialises global matrix
S = zeros(No,1);    % Initialises solution vector

switch mode % Switches between original source term and modified
    
    case 1 % Original source term (section 1 - 2a)
        
        %Loops for all element in matrix
        for eID = 1:Ne
            
            % Creates local element matrix for diffusion operator
            elematLap = LaplaceElemMatrix(D,eID,mesh);
            
            % Creates local element matrix for linear reaction operator
            elematLin = LinOpElemMatrix(L,eID,mesh);
            
            % Calculates overall local element matrix
            elementLoc = elematLap - elematLin;
            
            No = eID+1;
            
            % Sums local element matrix into global matrix accross elements
            GM(eID:No,eID:No) = GM(eID:No,eID:No) + elementLoc;
            
            % Returns the Jacobian from the data structure
            J = mesh.elem(eID).J;
            
            % Adds source term element vector to solution vector
            S(eID:No,1) = S(eID:No,1) + (f*J);
            
        end
        
    case 2 % Modified source term (question 2b)
        
        %Loops for all element in matrix
        for eID = 1:Ne
            
            % Creates local element matrix for diffusion operator
            elematLap = LaplaceElemMatrix(D,eID,mesh);
            
            % Creates local element matrix for linear reaction operator
            elematLin = LinOpElemMatrix(L,eID,mesh);
            
            % Calculates overall local element matrix
            elementLoc = elematLap - elematLin;
            
            No = eID+1;
            
            % Sums local element matrix into global matrix accross elements
            GM(eID:No,eID:No) = GM(eID:No,eID:No) + elementLoc;
            
            % Returns the Jacobian from the data structure
            J = mesh.elem(eID).J;
            
            % Returns local element x values
            x = mesh.elem(eID).x;
            
            % Creates modified local element vector with modified source
            % term
            Int0 = J*f*(1 + (8/3)*x(1)+(4/3)*x(2));
            Int1 = J*f*(1 + (4/3)*x(1)+(8/3)*x(2));
            
            % Adds source term element vector to solution vector
            S(eID:No,1) = S(eID:No,1) + ([Int0; Int1]);
            
        end
        
end

end
