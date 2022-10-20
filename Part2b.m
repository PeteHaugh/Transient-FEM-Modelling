function [] = Part2b(G)
% This function runs section 2b of the assignment

disp('Simultation in progress')

% Set Inputs
Ne = 60;
No = (Ne * 2) + 1;
% X values for stability chart
X.Val = 1/6;
X.MN = round(X.Val  * No);

% Time data structure
Time.dt = 0.05;
Time.N = 1000;
Time.theta = 1;
T = Time.N * Time.dt;

% Set lambda and solution to 0
MC.L = 0;
MC.f = 0;

% Set boundary condition
BC = [393.15, 310.15];

%Initial conditions
IC = [310.15, 310.15];

% Set domin
MC.xmin = 0;
MC.xmax = 0.01;

% Solve transient problem
[C,TS,~] = TransientSolverPart2(MC,Time,BC,IC,Ne,G);

% Find the element number where a burn occurs
TBN = find(C(X.MN,:) > 317.15);

% Find temperature value of burn (TB)
TB = TS(TBN(2));
x = round(TB/Time.dt);
TE = C(X.MN,x:end)';
integrand = (2 * 10 ^ 98) * exp( - (12017./(TE-273.15)));

%format short
Gamma = Time.dt * trapz(integrand);

if Gamma >= 1
    h = ['Second degree burn will occur',' (Gamma = ',num2str(round(Gamma,2,'significant')),')'];
else
    h = ['Second degree burn will notoccur','(Gamma = ',num2str(round(Gamma,2,'significant')),')'];
end
disp(h)
end




