function [F,M,K,mesh]= GlobalMatrix(MC,Ne)
% Global Matrix Creator for part 2
% This function creates the global matrix for part 1 of the assignment

No = 2 * Ne + 1; 

M = zeros(No);      % Initialises mass matrix
K = zeros(No);      % Initialises stiffness matrix
F = zeros(No, 1);

[mesh] = OneDimLinearMeshGen(MC.xmin,MC.xmax,Ne);

% Run loop to allocate different material parameters for each skin layer
for eID = 1:Ne
    
    [LocalM] = MassElemMatrix(eID,mesh); % Creates local mass matrix
    
    % Global matrix construction logic 
    a = eID * 2 - 1;
    b = a + 2;
    
    [elematLap] = LaplaceElemMatrix(MC.D,eID,mesh); % Creates local diffusion coefficient
 
    [elematLin] = LinElemMatrix(MC.L,mesh); % Creates local linear reaction matrix
    
    [elematSol] = SolElemMatrix(MC.f,eID,mesh); % Creates local solution vector
    
    M(a:b,a:b) = M(a:b,a:b) + LocalM; % Constructs mass matrix
    
    elematLoc = elematLap - elematLin;
    
    K(a:b,a:b) = K(a:b,a:b) + elematLoc; % Constructs stiffness matrix
    
    F(a:b,1) = F(a:b,1) + elematSol; % Constructs solution vector
    
end

end