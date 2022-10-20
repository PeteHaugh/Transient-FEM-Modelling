function [] = Part2c(G)
% This function runs section 2c of the assignment

disp('Simultation in progress')

% Set Inputs
Ne = 60;
No = (Ne * 2) + 1;

% X values for stability chart
X.Val = 1/6;
X.MN = round(X.Val * No);

% Time data structure
Time.dt = 0.05;
Time.N = 1000;
Time.theta = 1;

% Set lambda and solution to 0
MC.L = 0;
MC.f = 0;

% Set boundary condition
BC = [393.15, 310.15];

%Initial conditions
IC = [310.15, 310.15];

% Set material coefficient
MC.xmin = 0;
MC.xmax = 0.01;

% Arbritrary value to start while loop
Gamma = 2;

while Gamma > 1
    
    [C,TS,~] = TransientSolverPart2(MC,Time,BC,IC,Ne,G);
    
    TBN = find(C(X.MN,:) > 317.15);
    TB = TS(TBN(2));
    x = round(TB/Time.dt);
    
    TE = C(X.MN,x:end)';
    
    integrand = (2 * 10 ^ 98) * exp( - (12017./(TE-273.15)));
    
    Gamma = Time.dt * trapz(integrand);
    
    if Gamma > 1 * 10^10
        
        BC(1) = BC(1) - 10;
        
    elseif Gamma < 1 * 10^10 && Gamma > 1 * 10^5
        
        BC(1) = BC(1) - 1;
        
    elseif Gamma < 1 * 10^5
        
        BC(1) = BC(1) - 0.5;
    end
    
    disp(['Gamma value at ', num2str(BC(1)),'K outer BC is ',num2str(round(Gamma,3,'significant'))])
    
end

Reduc = 393.15 - BC(1);
disp(['Final boundary condition = ',num2str(round(BC(1),1)),'K'])
disp(['Minimum temperature reduction of ', num2str(round(Reduc,1)), 'K required'])

end