function [] = TransientStability

% Make into function with T variable

Ne = 17;

No = 2 * Ne + 1;        % Converts No. of elements to No. of nodes

T = 1;
theta = 0.5;
dt = 0.1;
N = T/dt;
xmin = 0;
xmax = 1;
D = 1;
L = 0;
f=0;

Ccurr = zeros(No,1);
S = zeros(No,1);

%Initial conditions
figure;
for t = 1:N
    
    [M,K,mesh] = GlobalMatrix(D,L,xmin,xmax,Ne);
    
    GMcurr = M + (theta * dt * K);
    
    GMprev = M - (1 - theta) * dt * K;
    
    GV = GMprev * Ccurr;
    
    for eID = 1:Ne-1
        
        S(eID) = dt*(theta-1)*S(eID);   
        S(eID+1) = dt*theta*S(eID+1);
        
        NBC(eID) = dt*(theta-1)*NBC(eID);   
        NBC(eID+1) = dt*theta*NBC(eID+1);
        
    end
    
    GV = GV + S + NBC;
    
    GV(1) = 0;
    GV(end) = 1; % TEMPORARY BC FIX
    
    GMcurr(1,:) = zeros(1,No);
    GMcurr(No,:) = zeros(1,No);
    GMcurr(1,1) = 1;
    GMcurr(No,No) = 1;
    
    Cnext = GMcurr\GV;
    
    Ccurr = Cnext;
    
    %     if(mod(t,5)==0)
    %         x = mesh.nvec;
    %         CVec = Ccurr(1:2:end);
    %         plot(x,CVec);
    %
    %     end
    
    % When x = 0.8 plot results
    % Store as vector
    
    % 80% along vector and round up or interpolate
    
    CVec = Ccurr(1:2:end);
    
    
end

x = mesh.nvec


figure;
hold on

title('FIXED THIS FUCKING FUCKED FUCK')
legend('FUCK THIS LINE')
