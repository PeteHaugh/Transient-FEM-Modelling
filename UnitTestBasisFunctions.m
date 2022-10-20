%% Test 1:test that psi(1) is the correct value
tol = 1e-14;

gq = CreateGQScheme(3); 

[psi] = EvalBasis(1,gq.xipts);

% Test that the correct psi value is calculated
assert((abs(psi(1) - (-sqrt(3/5)*(-sqrt(3/5)-1)/2))) <= tol)
assert((abs(psi(2) - (0*(0-1)/2))) <= tol)
assert((abs(psi(3) - (sqrt(3/5)*(sqrt(3/5)-1)/2))) <= tol)

%% Test 2:test that psi(2) is the correct value
tol = 1e-14;

gq = CreateGQScheme(3); 

[psi] = EvalBasis(2,gq.xipts);

% Test that the correct psi value is calculated
assert((abs(psi(1) - (1 - 3/5))) <= tol)
assert((abs(psi(2) - (1 - 0))) <= tol)
assert((abs(psi(3) - (1 - 3/5))) <= tol)

%% Test 3:test that psi(3) is the correct value
tol = 1e-14;

gq = CreateGQScheme(3); 

[psi] = EvalBasis(3,gq.xipts);

% Test that the correct psi value is calculated
assert((abs(psi(1) - (-sqrt(3/5)*(-sqrt(3/5)+1)/2))) <= tol)
assert((abs(psi(2) - (0*(0+1)/2))) <= tol)
assert((abs(psi(3) - (sqrt(3/5)*(sqrt(3/5)+1)/2))) <= tol)

%% Test 4:test that dpsidxi(1) is the correct value
tol = 1e-14;

gq = CreateGQScheme(3); 

dpsidxi = EvalBasisGrad(1,gq.xipts);

% Test that the correct dpsidxi value is calculated
assert((abs(dpsidxi(1) - (-sqrt(3/5) - 0.5))) <= tol)
assert((abs(dpsidxi(2) - (0 - 0.5))) <= tol)
assert((abs(dpsidxi(3) - (sqrt(3/5) - 0.5))) <= tol)


%% Test 5:test that dpsidxi(2) is the correct value
tol = 1e-14;

gq = CreateGQScheme(3); 

dpsidxi = EvalBasisGrad(2,gq.xipts);

% Test that the correct dpsidxi value is calculated
assert((abs(dpsidxi(1) - (-2 * -sqrt(3/5)))) <= tol)
assert((abs(dpsidxi(2) - (-2 * 0))) <= tol)
assert((abs(dpsidxi(3) - (-2 * sqrt(3/5)))) <= tol)

%% Test 6:test that dpsidxi(3) is the correct value
tol = 1e-14;

gq = CreateGQScheme(3); 

dpsidxi = EvalBasisGrad(3,gq.xipts);

% Test that the correct dpsidxi value is calculated
assert((abs(dpsidxi(1) - (-sqrt(3/5) + 0.5))) <= tol)
assert((abs(dpsidxi(2) - (0 + 0.5))) <= tol)
assert((abs(dpsidxi(3) - (sqrt(3/5) + 0.5))) <= tol)

%% Test 6:test that GQ scheme is correct

tol = 1e-14;

gq = CreateGQScheme(3); 

% Check that correct values are obtained
assert((abs(gq.gsw(1) - (5/9)) <= tol) && (abs(gq.gsw(2) - (8/9)) <= tol) && (abs(gq.gsw(3) - (5/9)) <= tol))

assert((abs(gq.xipts(1) + sqrt(3/5)) <= tol) && (abs(gq.xipts(2) - 0) <= tol) && (abs(gq.xipts(3) - sqrt(3/5)) <= tol))


