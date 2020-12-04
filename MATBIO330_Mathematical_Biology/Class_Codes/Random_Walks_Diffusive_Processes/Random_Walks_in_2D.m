%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes Random Walks in 2D to compute root-mean-squared
%           distance from starting point.
%
%
% Author: Nick Battista
% Institution: TCNJ
% Created: April 8, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Random_Walks_in_2D()

M = 10000; % # of Random Walkers

N = 50;  % # of steps for each walker

ds = 0.1; % size of step (step-size)

x = 0;    % initial x-Position
y = 0;    % initial y-Position

% Perform a Random Walk for Each Walker
for i=1:M
   
    % Do Random Walk for i-th random walker
    [xFinal,yFinal,rSqr] = do_Random_Walk(N,ds,x,y);
    
    % Store x,y-Final positions
    xFinal_Vec(i) = xFinal;
    yFinal_Vec(i) = yFinal;
    
    % Store squared displacement distance from starting point (must be +)
    rSqr_Vec(i) = rSqr;
    
end


% Compute Avg. of Squared-Displacement Distance from Starting Pt.
RMS = sqrt( mean( rSqr_Vec ) );


% Print Information to Screen for RMS (Room Mean Squared-Displacement) From Starting Pt.
fprintf('Avg. Squared-Displacement: %2.4f\n',RMS);
fprintf('Theory Says Avg. Displacement = %2.4f\n',sqrt(N)*ds);
fprintf('Error: %2.4f\n\n\n',abs(RMS - sqrt(N)*ds));

% Plot ending point (x,y) for each Random Walk
figure(2)
plot(xFinal_Vec,yFinal_Vec,'b.','MarkerSize',10); hold on;
plot(0,0,'r.','MarkerSize',50); hold on;
xlabel('x');
ylabel('y');
title('2D Random Walks Final Position');


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

function [xFinal,yFinal,rSqr] = do_Random_Walk(N,ds,x,y)


% Perform the Random Walk
for i=1:N
   
    % Get random angle between 0 and 2*pi
    ang = 2*pi*rand(1); 
    
    % Move in x-direction (based on trig relations)
    x = x + ds*cos(ang);
    
    % Move in y-direction (based on trig relations)
    y = y + ds*sin(ang);
        
end

% Define final positions of (x,y) for output
xFinal = x;
yFinal = y;

% Save r^2 value
rSqr = x^2 + y^2;




