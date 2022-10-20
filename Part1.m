function [] = Part1(theta)
% This function runs part one of the assignment

% Set Inputs
Ne = 10;
No = (2 * Ne) + 1;

% Time data structure
Time.dt = 0.01;
Time.N = 100;
Time.theta = theta;

% Set material coefficient data structure
MC.D = 1;
MC.L = 0;
MC.f = 0;
MC.xmin = 0;
MC.xmax = 1;

% Set boundary condition
BC = [0, 1];

% Initial conditions
IC = [0, 0];

% Runs TransientSolver.m and solves for given coefficients and conditions
[C, TS, AnaC, ~, ~, X, Error] = TransientSolver(MC,Time,BC,IC,Ne);

% Creates a a vector containing evenly disptributed x values depending on
% the number of points
xvec = zeros(No, 1);

for x = 1 : No
    
    dx = MC.xmax / No;
    xvec(x) = (x - 1) * dx;
    
end

% Plots the solution against spatial distance in the domain
hfig1 = figure(1);
hold on
plot(xvec, C(:,(round(0.05/Time.dt)+1)),'-o')
plot(xvec, C(:,(round(0.1/Time.dt)+1)),'-o')
plot(xvec, C(:,(round(0.3/Time.dt)+1)),'-o')
plot(xvec, C(:,(round(1/Time.dt)+1)),'-o')
grid on
title('Plot of Solution C Against Distance x','FontSize', 24)
xlabel('Distance x','FontSize', 18)
ylabel('Solution c','FontSize', 18)
legend('t = 0.05s', 't = 0.10s', 't = 0.30s','t = 1.00s','Location','best','FontSize', 18)
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',18);
saveas(hfig1,'Part1a.png')

% Plots the solution against time at a distance of x = 0.8
hfig2 = figure(2);
hold on
plot(TS, C(X.MN,:),'-o')
plot(TS, AnaC)
grid on
axis([0 1 0 1])
title('Solution C Against Time t at Tistance x = 0.8','FontSize', 24) 
xlabel('Time t (s)','FontSize', 18)
ylabel('Solution c','FontSize', 18)
legend('Numerical','Analytic','FontSize', 18)
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',18);
saveas(hfig2,'Part1b.png')

% Plots the error between numerical and analytical value at a distance of x
% = 0.8m
hfig3 = figure(3);
hold on
plot(TS, Error)
grid on
title('Error Against Time t at Distance x = 0.8','FontSize', 24) 
xlabel('Time t (s)','FontSize', 18)
ylabel('Error','FontSize', 18)
txt = ['Max error = ',num2str(max(Error))];
dim = [.3 .5 .3 .3];
annotation('textbox',dim,'String',txt,'FontSize',24,'FitBoxToText','on');
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',18);
saveas(hfig3,'Error.png')






