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

function Random_Walks_in_2D_Lattice()

M = 10000; % # of Random Walkers

N = 300;  % # of steps for each walker

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
RMS = sqrt( rSqr_Avg );


% Print Information to Screen for RMS (Room Mean Squared-Displacement) From Starting Pt.
fprintf('Avg. Squared-Displacement: %2.4f\n',RMS);
fprintf('Theory Says Avg. Displacement = %2.4f\n',sqrt(N)*ds);
fprintf('Error: %2.4f\n\n\n',abs( RMS - sqrt(N)*ds));

% Plot ending point (x,y) for each Random Walk
plot(xDist_Vec,yDist_Vec,'b.','MarkerSize',10); hold on;
plot(0,0,'r.','MarkerSize',60); hold on;
xlabel('x');
ylabel('y');
title('2D Random Walks Final Position');
%axis square;
axis([-5 5 -5 5]);
set(gca,'FontSize',18);
grid on;

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




