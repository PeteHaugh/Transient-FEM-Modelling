%% Test 1: test symmetry of the global matrix
% Test that this global matrix is symmetric
tol = 1e-14;

% Material coefficients are defined
MC.xmin = 0;       % Sets Domain
MC.xmax = 1;
Ne = 10;           % Sets number of elements
MC.D = 1;          % Diffusion coefficient
MC.L = 1;          % Linear coefficient
MC.f = 0;

[~,M,K,~] = GlobalMatrix(MC,Ne);


assert(abs(K(1,2) - K(2,1)) <= tol)
assert(abs(M(1,2) - M(2,1)) <= tol)
assert(abs(K(2,3) - K(3,2)) <= tol)
assert(abs(M(2,3) - M(3,2)) <= tol)

%% Test 2: test global matrix values are correctly summed
% % Tests that the values of the global matrix are being summed correctly

tol = 1e-14;

% Material coefficients are defined
MC.xmin = 0;       % Sets Domain
MC.xmax = 1;
Ne = 10;           % Sets number of elements
MC.D = 1;          % Diffusion coefficient
MC.L = 1;          % Linear coefficient
MC.f = 0;

[~,M,K,~] = GlobalMatrix(MC,Ne);

diff = M(1,1) - 0.5 * M(3,3);
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)

diff = K(1,1) - 0.5 * K(3,3);
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)

%% Test 3: test solution vector values are correctly summed
% % Tests that the values of the solution vector are being summed correctly
tol = 1e-14;

% Material coefficients are defined
MC.xmin = 0;       % Sets Domain
MC.xmax = 1;
Ne = 10;           % Sets number of elements
MC.D = 1;          % Diffusion coefficient
MC.L = 1;          % Linear coefficient
MC.f = 0;

[F,~,~,~] = GlobalMatrix(MC,Ne);

diff = F(1) - 0.5 * F(3);
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)


%% Test 4: test that Dirichlet boundary conditions are correctly applied
% % Tests that Dirichlet boundary conditions are correctly applied to the
% % top and bottom rows of the global matrix

tol = 1e-5;

% Material coefficients are defined
MC.xmin = 0;       % Sets Domain
MC.xmax = 1;
Ne = 10;           % Sets number of elements
No = 2 * Ne + 1;
MC.D = 1;          % Diffusion coefficient
MC.L = 1;          % Linear coefficient
MC.f = 0;

% Time data structure
Time.dt = 0.01;
Time.N = 100;
Time.theta = 1;

t = 0.5;

MN = round((t/Time.dt) + 1);

% Set boundary condition
BC = [0, 1];

C = zeros(No, Time.N + 1);

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

assert(abs(GM(1,1) - 1) <= tol)
assert(GM(No,No) - 1 <= tol)

for a = 2:No
    assert(GM(1,a) <= tol)
end

for a = 1:No-1
    assert(GM(No,a) <= tol)
end
