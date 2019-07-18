%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes error btwn simulation and theory for different numbers
%           of random walkers in 2D. As # of RWs goes up, error goes down.
%
% Author: Nick Battista
% Institution: TCNJ
% Created: April 8, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compute_Random_Walks_2D_Lattice_Convergence()

% Vector of all #'s of Random Walker Trials to Do
MVec = [10:10:90 100:100:900 1e3:1e3:9e3 1e4:1e4:9e4 1e5:1e5:5e5];

for i=1:length(MVec)
    
    % Define # of Random Walkers for Trial
    M = MVec(i);
    
    % Print Simulation Case Info to Screen
    fprintf('Simulation Case w/ M=%d\n',M);
    
    % Call Random Walk Function that Returns Error btwn theory and simulation for RMS
    err = Random_Walks_in_2D_Lattice(M);
    
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

function errRMS = Random_Walks_in_2D_Lattice(M)

%M = 10000; % # of Random Walkers

N = 25;  % # of steps for each walker

ds = 0.1; % size of step (step-size)

x = 0;    % initial x-Position
y = 0;    % initial y-Position

% Perform a Random Walk for Each Walker
for i=1:M
   
    % Do Random Walk for i-th random walker
    [xDist,yDist,rSqr] = do_Random_Walk(N,ds,x,y);
    
    % Store displacement distance from starting point (can be + or -)
    xDist_Vec(i) = xDist;
    yDist_Vec(i) = yDist;

    % Store squared displacement distance from starting point (must be +)
    rSqr_Vec(i) = rSqr;
    
end


% Compute Avg. of Squared-Displacement Distance from Starting Pt.
rSqr_Avg = mean(rSqr_Vec);

% Compute root-mean-squared displacement (take sqrt of avg. r-Squared Displacement Vector)
RMS_Avg = sqrt( rSqr_Avg );

% Store RMS Error between average RMS from simulation and theory
errRMS = abs( RMS_Avg - sqrt(N)*ds);

% Print Information to Screen for RMS (Room Mean Squared-Displacement) From Starting Pt.
% fprintf('Avg. Squared-Displacement: %2.4f\n',RMS_Avg);
% fprintf('Theory Says Avg. Displacement = %2.4f\n',sqrt(N)*ds);
% fprintf('Error: %2.4f\n\n\n',errRMS);

% Plot ending point (x,y) for each Random Walk
% plot(xDist_Vec,yDist_Vec,'b.','MarkerSize',10); hold on;
% plot(0,0,'r.','MarkerSize',50); hold on;
% xlabel('x');
% ylabel('y');
% title('2D Random Walks Final Position');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: perform Random Walk starting at "x,y" and ending after "N" steps
%
% Inputs: x,y  <-- starting point for x,y
%         N    <-- length of Random Walk
%         ds   <-- length of a step 
%
% Outputs:
%         xDist: displacement distance from starting point (can be + or -)
%         xSqr: squared displacement distance from starting point (must be +)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xDist,yDist,rSqr] = do_Random_Walk(N,ds,x,y)


% Perform the Random Walk
for i=1:N
   
    coin = rand(1); 
    
    if (coin <= 0.25)
       
        x = x - ds;
        
    elseif (coin <= 0.5)
        
        x = x + ds;
        
    elseif (coin <=0.75)
        
        y = y - ds;
        
    else
        
        y = y + ds;
    end
        
end

% Store final positions of (x,y) for Random Walker 
xDist = x;
yDist = y;

% Store squared displacement distance from origin  
rSqr = x^2 + y^2;




