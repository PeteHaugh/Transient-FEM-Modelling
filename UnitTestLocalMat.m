%% Test 1: test 2 different elements of the same size produce same matrix
% % Test that for two elements of an equispaced mesh, as described in the
% % lectures, the element matrices calculated are the same
tol = 1e-14;
D = 5; %diffusion coefficient
eID=1; %element ID
mesh = OneDimLinearMeshGen(0,1,10);

elemat1 = LaplaceElemMatrix(D,eID,mesh);

eID=2; %element ID

elemat2 = LaplaceElemMatrix(D,eID,mesh);

diff = elemat1 - elemat2;
diffnorm = sum(sum(diff.*diff));
assert(abs(diffnorm) <= tol)

%% Test 2: test that one matrix is evaluted correctly
% % Test that element 1 of the three element mesh problem described in the lectures
% % the element matrix is evaluated correctly
tol = 1e-14;
D = 1; %diffusion coefficient
eID=1; %element ID
mesh = OneDimLinearMeshGen(0,1,3);

elemat1 = LaplaceElemMatrix(D,eID,mesh);

elemat2 = [ 7 -8 1; -8 16 -8;1 -8 7];
diff = elemat1 - elemat2; %calculate the difference between the two matrices
diffnorm = sum(sum(diff.*diff)); %calculates the total squared error between the matrices
assert(abs(diffnorm) <= tol)

%% Test 3: test Symmetry of the Matrix
% Test to check matrix symmetry for the diffusion and linear reaction operators
tol = 1e-14;
D = 1; %diffusion coefficient
eID=1; %element ID
mesh = OneDimLinearMeshGen(0, 1, 10); 

elemat = LaplaceElemMatrix(D, eID, mesh); 
                      
assert((abs(elemat(1, 2) - elemat(2, 1)) <= tol))
assert((abs(elemat(1, 3) - elemat(3, 1)) <= tol)) 
assert((abs(elemat(2, 3) - elemat(3, 2)) <= tol))

elemat = LinElemMatrix(D, mesh); 
                      
assert((abs(elemat(1, 2) - elemat(2, 1)) <= tol))
assert((abs(elemat(1, 3) - elemat(3, 1)) <= tol)) 
assert((abs(elemat(2, 3) - elemat(3, 2)) <= tol))

