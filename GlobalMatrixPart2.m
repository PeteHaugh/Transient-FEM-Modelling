function [F,M,K,mesh]= GlobalMatrixPart2(MC,Ne,G)
% Global Matrix Creator for part 2
% This function creates the global matrix for part 2 of the assignment

%Set number of nodes
No = 2 * Ne + 1; 

M = zeros(No);      % Initialises mass matrix
K = zeros(No);      % Initialises stiffness matrix
F = zeros(No, 1);   % Initialises solution matrix

%Initialise mesh
[mesh] = OneDimLinearMeshGen(MC.xmin,MC.xmax,Ne);

% Set skin layer coordinates in mesh
E = Ne/6;
D = Ne/2;
B = Ne;

% Run loop to allocate different material parameters for each skin layer
for eID = 1:Ne
    
    if eID <= (E)
        k = 25;
        rho = 1200;
        c = 3300;
        
        MC.D = k / (rho * c);
        MC.L = 0;
        MC.f = 0;
        
    elseif ((E + 1) <= eID && (eID <= D))
        k = 40;
        rho = 1200;
        c = 3300;
        rhob = 1060;
        cb = 3770;
        Tb = 310.15;
        
        MC.D = k / (rho * c);
        MC.L = -(G * rhob *cb / (rho * c));
        MC.f = Tb * G * rhob *cb / (rho * c);
        
    elseif ((D + 1) <= eID && (eID <= B))
        k = 20;
        rho = 1200;
        c = 3300;
        rhob = 1060;
        cb = 3770;
        Tb = 310.15;
        
        MC.D = k / (rho * c);
        MC.L = -(G * rhob *cb / (rho * c));
        MC.f = Tb * G * rhob *cb / (rho * c);
    end
    
    [LocalM] = MassElemMatrix(eID,mesh); % Creates local mass matrix
    
    % Global matrix construction logic 
    a = eID * 2 - 1;
    b = a + 2;
    
    [elematLap] = LaplaceElemMatrix(MC.D,eID,mesh); % Creates local diffusion coefficient matrix
 
    [elematLin] = LinElemMatrix(MC.L,mesh); % Creates local linear reaction matrix
    
    [elematSol] = SolElemMatrix(MC.f,eID,mesh); % Creates local solution vector
    
    M(a:b,a:b) = M(a:b,a:b) + LocalM; % Constructs mass matrix
    
    elematLoc = elematLap - elematLin;
    
    K(a:b,a:b) = K(a:b,a:b) + elematLoc; % Constructs stiffness matrix
    
    F(a:b,1) = F(a:b,1) + elematSol; % Constructs solution vector
    
end

end