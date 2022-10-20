% Assignment 1 Master Script
% % This script allows selectin of the required report section

clear
close

qlist = {'1', '2a','2b','2c','2d'}; 
[question] = listdlg('PromptString','Please enter quesiton number:','liststring',qlist); %Handle response

switch question
    case 1
        % Runs section 1: Verifies software against analytical solution and runs unit testing
        
        qlist = {'Crank-Nicolson', 'Forward Euler','Backward Euler'}; 
        [TSM] = listdlg('PromptString','Please choose timestepping method','liststring',qlist); %Handle response
        
        switch TSM
            case 1 % Crank-Nicolson
                theta = 0.5;
            case 2 % Forward Euler
                theta = 0;
            case 3 % Backward Euler
                theta = 1;   
        end
        
        Part1(theta)
        
        L2NormError
        
        rt1 = table(runtests('UnitTestBasisFunctions'));
        disp(rt1);
        rt2 = table(runtests('UnitTestLocalMat'));
        disp(rt2);
        rt3 = table(runtests('UnitTestGlobalMatrix'));
        disp(rt3);
         
    case 2
        % Runs section 2a showing spatial temperture distribution in the tissue between t=0 and t=50s
        G = 0;
        Part2a(G)
        
    case 3
        % Runs Section 2b showing whether a second degree burn will occur in the epidermis layer after t=50s 
        G = 0;
        Part2b(G)
        
    case 4
        % Runs section 2c which determines the minimum temperature reduction required by protective clothing so that no second degree burn occurs
        G = 0;
        Part2c(G)
        
    case 5
        % Reruns section  2a and 2c but including blood flow related reaction & source terms.
        G = 0.0375;
        Part2d(G)
        Part2c(G)
        
end





