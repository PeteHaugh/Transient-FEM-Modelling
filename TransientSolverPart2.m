function [C, TS, mesh] = TransientSolverPart2(MC,Time,BC,IC,Ne,G)
% Transient problem solver
% This function solves transient problems over a given time period 

% Calculate simulation time from time data structure
T = Time.dt * Time.N;

% Calculate number of nodes
No = 2 * Ne + 1;  

C = zeros(No, Time.N + 1);       % Initialise solution vector
TS = zeros((Time.N + 1), 1);     % Initialise timestep vector

% Set IC in solution vector
C(:, 1) = IC(2);
C(1, 1) = IC(1);

% Initialise global matrix (GM) and vector (GV)
GM = zeros(No, No);
GV = zeros(No, 1); 

for t = 0:Time.dt:T  
    
    % Calculate and store the matrix number of time step
    MN = round((t/Time.dt) + 1);
    
    % Create source (F), Mass (M) and stiffness (K) matrices and initialise
    % mesh
    [F,M,K,mesh] = GlobalMatrixPart2(MC,Ne,G);
    
    % Calculate GM values
    GM = M + (Time.theta * Time.dt * K);
    
    % Calculate GV values
    GV = (M - ((1 - Time.theta) * Time.dt * K)) * C(:,MN);
    
    % Calculate previous matrix and multiply by previous solution, then
    % store in GV
    GV = GV + (Time.dt * ((Time.theta * F) + ((1 - Time.theta) * F)));
    
    % Set Dirichlet boundary conditions
    [GM,GV] = Dirichlet(No,GM,GV,BC);
    
    % Solve final matrix
    C(:,MN + 1) = GM\GV;
    
    % Record the time step value 
    TS(MN) = t;
   
end

% Remove final entry
C = C(:,1:MN);

end