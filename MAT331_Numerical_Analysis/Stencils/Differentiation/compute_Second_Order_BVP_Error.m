%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nick Battista
% Instutition: TCNJ
% Course: MAT 331 (Numerical Analysis)
% Date: 3/4/18
%
% FUNCTION: compute the error when approximating the solution to the 
%           following BOUNDARY VALUE PROBLEM:
%           u'' = f
%           u(0) = 1; u(1) = 2;
%
%           NOTE: 1. f is constructed to already know the true solution. 
%                 2. u(x) = (1-x)*cos(x)+2xsin(0.5pi*x)
%                 3. f = 2sin(x) + (x-1)cos(x) + 2pi*cos( 0.5pi*x ) - 0.5pi^2*xsin(0.5pi*x)
%
% Inputs:  
%
% Returns:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compute_Second_Order_BVP_Error()

% How many boundary points to try
NVec = [5 10:10:90 100:100:900 1000:1000:5000];

% Initialize error vector / time storage vector
errVec = zeros(1,length(NVec));
time = errVec;

% Find error for each approximate solution
for i=1:length(NVec)
    
    % Compute error for approximation using Nvec(i) # of grid points
    tic
    errVec(i) = Second_Order_BVP( NVec(i) );
    
    time(i) = toc;
end

% Calculate Approximate Order of Convergence!
calculate_Order_Of_Convergence(errVec,NVec);

% PLOTS ERROR VS. GRID RESOLUTION (N)
figure(1) 
ms = 30; lw = 5;
loglog(NVec,errVec,'.','MarkerSize',ms,'LineWidth',lw);
xlabel('N (grid resolution)');
ylabel('Max. Abs. Error');
set(gca,'FontSize',18)

% PLOTS HOW LONG IT TAKES TO RUN FOR EACH N
figure(2) 
loglog(NVec,time,'.','MarkerSize',ms,'LineWidth',lw);
xlabel('N (grid resolution)');
ylabel('Time');
set(gca,'FontSize',18)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: calculates approximate order of convergence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calculate_Order_Of_Convergence(errVec,NVec)

% CALCULATE APPROXIMATE SLOPE:
slope = ( log(errVec(10)) - log(errVec(20)) ) / ( log(NVec(10)) - log(NVec(20)) );
fprintf('\n\nOrder of Convergence is approximately: %d\n\n',slope);
