%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Solves Discrete Dynamical Systems in Population Ecology
%
% Author: Nick Battista
% Institution: TCNJ
% Created: March 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Predator_Prey()

%
% Clears any previous plots that are open in MATLAB
clf;

%
% Time Information / Initialization
%
TFinal = 150;          % Simulation runs until time = TFinal
dt = 0.0025;           % Time-step
t = 0;                 % Initialize Time to 0.
n = 0;                 % Initialize storage counter to 0.


%
% Initial Values
%
x0 = 100;   % Initial Population for Prey, X
y0 = 2;    % Initial Population for Predator, Y


%
% Parameter Values
%
k = 0.75;    % growth rate
C = 250;     % carrying capacity
b1 = 0.075; % death parameter for prey from predator interactions
b2 = 0.05;   % growth parameter for predator from prey interactions
d = 0.5;       % death rate parameter for predator


%
% Saves Initial Values into Storage Vectors
%
X(1) = x0;
Y(1) = y0;
TimeVec(1) = t;

%
% While-loop that iteratively solves the discrete dynamical system
%
while t<TFinal
   
    % Iterate storage counter and time, t
    n = n + 1;
    t = t + dt;
    
    % Solve Prey ODE w/ Euler Method
    X(n+1) = X(n) + dt * ( k*X(n)*( 1 - X(n)/C ) - b1*X(n)*Y(n) ); 
    
    % Solve Predator ODE w/ Euler Method
    Y(n+1) = Y(n) + dt * ( -d*Y(n) + b2*X(n)*Y(n) );
    
    % Next time in time storage vector
    TimeVec(n+1) = t;
    
end



%
% Plot Attributes
%
lw = 4;  % LineWidth (how thick the lines should be)
ms = 25; % MarkerSize (how big the plot points should be)
fs = 18; % FontSize (how big the font should be for labels)

%
% PLOT 1: Populations vs. Time
%
figure(1)
plot(TimeVec,X,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(TimeVec,Y,'r.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Time');
ylabel('Population');
leg = legend('Prey','Predator');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);


%
% PLOT 2: Phase Plane Plot
%
figure(2)
plot(X,Y,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Prey Population');
ylabel('Predator Population');
set(gca,'FontSize',fs);
