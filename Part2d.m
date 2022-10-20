function [] = Part2d(G)
% This function runs part 2d of the assignment

% Set Inputs
Ne = 20;
No = (Ne * 2) + 1;

% X values for stability chart (Not used in this section)
X.Val = 1/6;
X.MN = round((X.Val * 2 * Ne) + 1);

% Time data structure
Time.dt = 0.01;
Time.N = 5000;
Time.theta = 1;

% Set lambda and solution to 0
MC.L = 0;
MC.f = 0;

% Set boundary condition
BC = [393.15, 310.15];

% Initial conditions
IC = [310.15, 310.15];

% Set material coefficient data structure
MC.xmin = 0;
MC.xmax = 0.01;

% Runs TransientSolver.m and solves for given coefficients and conditions
[C, TS, ~] = TransientSolverPart2(MC,Time,BC,IC,Ne,G);

% Re-runs section 2a) for comparison
G = 0;
[CG, TSG, ~] = TransientSolverPart2(MC,Time,BC,IC,Ne,G);

% Creates a a vector containing evenly disptributed x values depending on
% the number of points
xvec = zeros(No, 1);

for x = 1 : No
    
    dx = MC.xmax / No;
    xvec(x) = (x-1) * dx;
    
end

% Plots temperature throughout the skin layers
hfig = figure(1);
hold on
plot(xvec, C(:,(round(0.05/Time.dt) + 1)),'-o')
plot(xvec, CG(:,(round(0.05/Time.dt) + 1)),'-x')

plot(xvec, C(:,(round(1/Time.dt) + 1)),'-o')
plot(xvec, CG(:,(round(1/Time.dt) + 1)),'-x')

plot(xvec, C(:,(round(50/Time.dt) + 1)),'-o')
plot(xvec, CG(:,(round(50/Time.dt) + 1)),'-x')

grid on
title('Plot of temperature T against distance x','FontSize', 32)
xlabel('Distance x (m)','FontSize', 24)
ylabel('Temperature T (K)','FontSize', 24)
legend('t = 0.05s w/ flow', 't = 0.05s w/o flow', 't = 1.00s w/ flow','t = 1.00s w/o flow','t = 50.0s w/ flow','t = 50.0s w/o flow','Location','best','FontSize', 24)
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',24);
saveas(hfig,'Part2d.png')

% Plots the solution against time at a distance of x = E
hfig2 = figure(2);
hold on
plot(TS, C(X.MN,:),'-')
plot(TSG, CG(X.MN,:),'-')
grid on
title('Temperature T Against Time t at Distance x = E','FontSize', 24) 
xlabel('Time t (s)','FontSize', 18)
ylabel('Temperature T (K)','FontSize', 18)
legend('x = E w/o flow', 'x = E w/ flow','Location','best','FontSize', 24)
txt = ['Max temperature with blood flow = ', num2str(max(C(X.MN,:)))];
dim = [.3 .5 .3 .3];
annotation('textbox',dim,'String',txt,'FontSize',24,'FitBoxToText','on');
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',18);
saveas(hfig2,'Part2d2.png')

end


