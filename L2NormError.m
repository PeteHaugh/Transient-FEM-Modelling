function L2NormError
% L2 Norm calulcator
% This function determines the L2 norm and plots the results

% Set material coefficient data structure
MC.D = 1;
MC.L = 0;
MC.f = 0;
MC.xmin = 0;
MC.xmax = 1;

% Time data structure
Time.dt = 0.01;
Time.N = 100;
Time.theta = 1;
MC.D = 1;

% Set boundary condition
BC = [0, 1];

%Initial conditions
IC = [0, 0];

% Initialise Gaussian quadrature scheme
N=3;
[gq] = CreateGQScheme(N);

% Initialise L2 and dx matrices
dx = zeros(10,1);
L2 = zeros(10,1);


% Loops for a range of element numbers
for Ne = 10:10:100
    
    [C, ~, ~, mesh, xvec, ~, ~] = TransientSolver(MC,Time,BC,IC,Ne);
    
    % Initialise L2Norm matrix
    L2Norm = zeros(Ne,1);
    L2Sol = zeros(Time.N);
    
    % Loop for chosen time period
    for t = 1 : Time.N
        
        % Loop for chosen time period
        for eID = 1 : Ne
            
            % Returns Jacobian
            J = mesh.elem(eID).J; 
            
            % Calculates local error
            [E] = ErrorElemMatrix(xvec, C, eID, t);
            
            % Sums local element matrix as per formula
            L2elemat = sum(E .* gq.gsw .* J);
            
            % Stores summed value in vector in correct position
            L2Norm(eID) = L2elemat;
            
        end
        
        % Solve for entire domain
        L2Sol(t) = sqrt(sum(L2Norm))/(Ne^3);
        
    end
    
    % Sets vector to correct length
    L2Sol = L2Sol(1:Time.N);
    
    % Sums L2 vector and scales according to the number of elements 
    L2(Ne/10) = sum(L2Sol);
    dx(Ne/10) = mesh.dx;
    
end

% Calculates natural log of L2 value
lnE = abs(log(L2));

% Calculates natural log of dx value
lndx = abs(log(dx));

% Calculates mean gradient
m = mean(gradient(lnE,lndx));

% Plots L2 norm graph and displays gradient
hfig4 = figure(4);
plot(lndx,lnE)
grid on
title('L2 Norm','FontSize', 24)
xlabel('Ln h','FontSize', 18)
ylabel('Ln E','FontSize', 18)
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',18);
txt = ['Gradient = ',num2str(m)];
dim = [.3 .5 .3 .3];
annotation('textbox',dim,'String',txt,'FontSize',24,'FitBoxToText','on');

saveas(hfig4,'L2Norm.png')

end
