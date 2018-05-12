%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nick Battista
% Instutition: TCNJ
% Course: MAT 331 (Numerical Analysis)
% Date: 3/22/18
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

function compute_Second_Order_BVP_Error_Compare()

% How many boundary points to try
NVec = [7 10:10:90 100:100:900 1000:1000:2000];

% Initialize error vector / time storage vector
errVec = zeros(1,length(NVec));
time = errVec;

% Find error for each approximate solution
for i=1:length(NVec)
    
    % Compute error for approximation using Nvec(i) # of grid points
    tic
    errVec3pt(i) = Second_Order_BVP( NVec(i) );
    time3pt(i) = toc;
    
    tic
    errVec5pt(i) = Second_Order_BVP_5Pt_Stencil( NVec(i) );
    time5pt(i) = toc;

    tic
    errVecMix(i) = Second_Order_BVP_Mixed_3_5Pt_Stencil( NVec(i) );
    timeMix(i) = toc;
    
end


% Calculate Approximate Order of Convergence!
calculate_Order_Of_Convergence(errVec3pt,NVec,'3 pt');
calculate_Order_Of_Convergence(errVec5pt,NVec,'5 pt');


% PLOTS ERROR VS. GRID RESOLUTION (N)
figure(1) 
ms = 48; lw = 4;
loglog(NVec,errVec3pt,'b.-','MarkerSize',ms,'LineWidth',lw); hold on;
loglog(NVec,errVec5pt,'r.-','MarkerSize',ms,'LineWidth',lw); hold on;
xlabel('N (grid resolution)');
ylabel('Max. Abs. Error');
legend('3pt Stencil','5pt Stencil');
set(gca,'FontSize',18)
set(legend,'FontSize',18);
pause();

% add Mixed pt. stencil info!
calculate_Order_Of_Convergence(errVecMix,NVec,'Mixed 3 & 5 pt');
loglog(NVec,errVecMix,'k.-','MarkerSize',ms,'LineWidth',lw); hold on;
legend('3pt Stencil','5pt Stencil','Mixed 3&5 Pt Stencil');


% PLOTS HOW LONG IT TAKES TO RUN FOR EACH N
%figure(2) 
%loglog(NVec,time3pt,'b.-','MarkerSize',ms,'LineWidth',lw); hold on;
%loglog(NVec,time5pt,'r.-','MarkerSize',ms,'LineWidth',lw); hold on;
%loglog(NVec,timeMix,'k.-','MarkerSize',ms,'LineWidth',lw); hold on;
%xlabel('N (grid resolution)');
%ylabel('Time');
%legend('3pt Stencil','5pt Stencil','Mixed 3&5pt Stencil);
%set(legend,'FontSize',18);
%set(gca,'FontSize',18)

fprintf('\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: calculates approximate order of convergence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calculate_Order_Of_Convergence(errVec,NVec,str)

% CALCULATE APPROXIMATE SLOPE:
if strcmp(str,'3pt')
    slope = ( log(errVec(10)) - log(errVec(20)) ) / ( log(NVec(10)) - log(NVec(20)) );
else
    slope = ( log(errVec(5)) - log(errVec(8)) ) / ( log(NVec(5)) - log(NVec(8)) );
end

strPrint = ['\n\nOrder of Convergence is approximately: %d for ' str ' stencil\n'];
fprintf(strPrint,slope);
