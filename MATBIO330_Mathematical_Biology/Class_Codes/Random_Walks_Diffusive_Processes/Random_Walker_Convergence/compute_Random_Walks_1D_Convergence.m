%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes error btwn simulation and theory for different numbers
%           of random walkers in 1D. As # of RWs goes up, error goes down.
%
% Author: Nick Battista
% Institution: TCNJ
% Created: April 8, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compute_Random_Walks_1D_Convergence()

% Vector of all #'s of Random Walker Trials to Do
MVec = [10:10:90 100:100:900 1e3:1e3:9e3 1e4:1e4:9e4 1e5:1e5:5e5];

for i=1:length(MVec)
    
    % Define # of Random Walkers for Trial
    M = MVec(i);
    
    % Print Simulation Case Info to Screen
    fprintf('Simulation Case w/ M=%d\n',M);
    
    % Call Random Walk Function that Returns Error btwn theory and simulation for RMS
    err = Random_Walks_in_1D(M);
    
    % Store RMS Error for Simulation with M random walkers
    errVec(i) = err;
    
end

figure(2)
lw=4;
ms=24;
loglog(MVec,errVec,'-.','MarkerSize',ms,'LineWidth',lw);
xlabel('# of Random Walkers');
ylabel('RMS Error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: perform Random Walk with M-Random Walkers to compute avg. 
%           error in RMS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errRMS = Random_Walks_in_1D(M)

%M = 1000; % # of Random Walkers

N = 20;  % # of steps for each walker

dx = 0.1; % size of step (step-size)

x = 0;    % initial position

% Perform a Random Walk for Each Walker
for i=1:M
   
    % Do Random Walk for i-th random walker
    [xDist,xSqr] = do_Random_Walk(N,dx,x);
    
    % Store displacement distance from starting point (can be + or -)
    xDist_Vec(i) = xDist;
    
    % Store squared displacement distance from starting point (must be +)
    xSqr_Vec(i) = xSqr;
    
end

% Compute Avg. of Displacement Distances from Starting Pt.
xDist_Avg = mean(xDist_Vec);

% Compute Avg. of Squared-Displacement Distance from Starting Pt.
RMS_Avg = sqrt(mean(xSqr_Vec));

% Store RMS Error between average RMS from simulation and theory
errRMS = abs( RMS_Avg - sqrt(N)*dx);

% Print Information to Screen for Avg. Displacement From Starting Pt.
%fprintf('\n\nAvg. Displacement: %2.4f\n',xDist_Avg);
%fprintf('Theory Says Avg. Displacement = 0\n');
%fprintf('Error: %2.4f\n\n\n',abs(xDist_Avg));

% Print Information to Screen for RMS (Room Mean Squared-Displacement) From Starting Pt.
%fprintf('Avg. Squared-Displacement: %2.4f\n',RMS_Avg);
%fprintf('Theory Says Avg. Displacement = %2.4f\n',sqrt(N)*dx);
%fprintf('Error: %2.4f\n\n\n',errRMS);


% Plot ending point (x,y) for each Random Walk
%plot(0,0,'r.','MarkerSize',50); hold on;
%plot(xDist_Vec,zeros(length(xDist_Vec)),'b.','MarkerSize',10); hold on;
%xlabel('x');
%ylabel('y');
%title('1D Random Walks Final Position');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: perform Random Walk starting at "x" and ending after "N" steps
%
% Inputs: x  <-- starting point
%         N  <-- length of Random Walk
%         dx <-- length of a step 
%
% Outputs:
%         xDist: displacement distance from starting point (can be + or -)
%         xSqr: squared displacement distance from starting point (must be +)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xDist,xSqr] = do_Random_Walk(N,dx,x)


% Perform the Random Walk
for i=1:N
   
    coin = rand(1); 
    
    if coin > 0.5
       
        x = x - dx;
        
    else
        
        x = x + dx;
    end
    
end

xDist = x;

xSqr = x^2;




