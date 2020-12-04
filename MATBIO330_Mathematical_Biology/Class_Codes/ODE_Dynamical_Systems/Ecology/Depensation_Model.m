%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Models a depensation differential equation based on the
%           Logistic Model in in Ecology
%
% Author: Nick Battista
% Institution: TCNJ
% Created: March 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Depensation_Model(x0)

%
% Clears any previous plots that are open in MATLAB
clf;

%
% Time Information / Initialization
%
TFinal = 30;           % Simulation runs until TFinal
dt = 0.0025;           % Time-step
t = 0;                 % Initialize Time to 0.
n = 0;                 % Initialize storage counter to 0.

%
% Initial Values
%
%x0 = 150;   % Initial Population
X(1) = x0;
tVec(1) = t;

%
% Parameter Values
%
k = 2;    % growth rate
C = 250;    % carrying capacity
r = 125;    % depensation constant

%
% While-loop that iteratively solves the differential equation
%
while t<TFinal
   
    % Iterate storage counter and time, t
    n = n + 1;
    t = t + dt;
    
    % Solve ODE w/ Euler Method
    X(n+1) = X(n) + dt* ( k*X(n)*( 1 - X(n)/C )*( X(n)/r - 1) ); 
    
    % Next time in time storage vector
    tVec(n+1) = t;
    
end

%
% Plot Attributes
%
lw = 4;  % LineWidth (how thick the lines should be)
ms = 25; % MarkerSize (how big the plot points should be)
fs = 18; % FontSize (how big the font should be for labels)

%
% PLOT 1: Population vs. Time
%
figure(1)
plot(tVec,X,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Time');
ylabel('Population');
set(gca,'FontSize',fs);
maxVal = 1.05*max(X);
axis([0 TFinal 0 maxVal]);

