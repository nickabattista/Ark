%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes Random Walks in 1D to compute root-mean-squared
%           distance from starting point.
%
%
% Author: Nick Battista
% Institution: TCNJ
% Created: April 8, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errRMS = Random_Walks_in_1D()

M = 1000; % # of Random Walkers

N = 50;  % # of steps for each walker

ds = 0.1; % size of step (step-size)

x = 0;    % initial position

% Perform a Random Walk for Each Walker
for i=1:M
   
    % Do Random Walk for i-th random walker
    [xDist,xSqr] = do_Random_Walk(N,ds,x);
    
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
errRMS = abs( RMS_Avg - sqrt(N)*ds);

% Print Information to Screen for Avg. Displacement From Starting Pt.
fprintf('\n\nAvg. Displacement: %2.4f\n',xDist_Avg);
fprintf('Theory Says Avg. Displacement = 0\n');
fprintf('Error: %2.4f\n\n\n',abs(xDist_Avg));

% Print Information to Screen for RMS (Room Mean Squared-Displacement) From Starting Pt.
fprintf('Avg. Squared-Displacement: %2.4f\n',RMS_Avg);
fprintf('Theory Says Avg. Displacement = %2.4f\n',sqrt(N)*ds);
fprintf('Error: %2.4f\n\n\n',errRMS);


% Plot ending point (x,y) for each Random Walk
plot(0,0,'r.','MarkerSize',50); hold on;
plot(xDist_Vec,zeros(length(xDist_Vec)),'b.','MarkerSize',10); hold on; % 'b.' <- blue, 'r.' <- red, 'g.' <- green, 'k.' <-black
xlabel('x');
ylabel('y');
title('1D Random Walks Final Position');
set(gca,'FontSize',18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: perform Random Walk starting at "x" and ending after "N" steps
%
% Inputs: x  <-- starting point
%         N  <-- length of Random Walk
%         ds <-- length of a step 
%
% Outputs:
%         xDist: displacement distance from starting point (can be + or -)
%         xSqr: squared displacement distance from starting point (must be +)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xDist,xSqr] = do_Random_Walk(N,ds,x)


% Perform the Random Walk
for i=1:N
   
    coin = rand(1); 
    
    if coin > 0.5
       
        x = x - ds;
        
    else
        
        x = x + ds;
    end
    
end

xDist = x;

xSqr = x^2;




