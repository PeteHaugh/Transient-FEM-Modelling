function [C, TS, AnaC, mesh, xvecL2, X, Error] = TransientSolver(MC,Time,BC,IC,Ne)
% Transient problem solver
% This function solves transient problems over a given time period 

% X values for stability chart
X.Val = 0.8;
X.MN = round((X.Val * 2 * Ne) + 1);

% Calculate number of nodes
T = Time.dt * Time.N;

No = 2 * Ne + 1;

C = zeros(No, Time.N + 1);       % Initialise solution vector
TS = zeros((Time.N + 1), 1);     % Initialise timestep vector
AnaC = zeros((Time.N + 1), 1);   % Initalise analytical solution vector
Error = zeros((Time.N + 1), 1);

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
    [F,M,K,mesh] = GlobalMatrix(MC,Ne);
    
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
    
    % Calculate analytical solution 
    AnaC(MN) = TransientAnalyticSoln(X.Val,t);
    
    % Calcualte error between numeric and analytical solution
    Error(MN) = abs(C(X.MN, MN) - AnaC(MN));
    
end

% Remove final entry
C = C(:,1:MN);

% Creates a a vector containing evenly disptributed x values depending on
% the number of points for running L2NormError.m function
xvecL2 = zeros(No, 1);

for x = 1 : (No + 1)
    
    mesh.dx = MC.xmax / ((2 * Ne));
    xvecL2(x) = (x-1) * mesh.dx;
    
end