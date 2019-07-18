%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Compares Discrete to Continuous Logistic Equation
%
% Author: Nick Battista
% Institution: TCNJ
% Created: March 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_Compare_Logistic(k)


%
% Clears any previous plots that are open in MATLAB
clf;

%
% Time Information / Initialization
%
TFinal = 100;          % Simulation runs until TFinal


%
% Initial Values
%
x0 = 25;


%
% Parameter Values
%
%k = 2.0;    % growth rate
C = 250;    % carrying capacity

%
% Call function to solve Discrete Dynamical System
%
[X_dis,tVec_Discrete] = please_Solve_Discrete_System(TFinal,k,C,x0);

%
% Call function to solve Continuous Dynamical System
%
[X_con,tVec_Continuous] = please_Solve_Continuous_System(TFinal,k,C,x0);

%
% Plot Attributes
%
lw = 3;  % LineWidth (how thick the lines should be)
ms = 30; % MarkerSize (how big the plot points should be)
fs = 18; % FontSize (how big the font should be for labels)

%
% PLOT 1: Populations vs. Time
%
figure(1)
plot(tVec_Discrete,X_dis,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(tVec_Continuous,X_con,'r-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Time');
ylabel('Population');
leg = legend('Discrete','Continuous');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: solves Logistic Differential Equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X,tVec] = please_Solve_Continuous_System(TFinal,k,C,x0)

%
% Initialize Time Information 
%
dt = 0.02;  % Time-step
t= 0;       % Initial Time

%
% Initialize Initial Population / Time Storage Vector / Counter for While Loop Indexing
%
X(1) = x0;
tVec(1) = 0;
n = 1;

%
% Solve ODE using Euler's Method
%
while t<TFinal
    
    % Update Time / Counter
    t = t + dt;
    n = n + 1;
    
    % Euler's Method to Solve for Solution
    X(n) = X(n-1) + dt*k*X(n-1)*( 1 - X(n-1) / C );
    
    % Update time-Vec
    tVec(n) = tVec(n-1) + dt;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: solves the Discrete Dynamical System for the Logistic Eqn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X,tVec] = please_Solve_Discrete_System(TFinal,k,C,x0)

%
% Initializing storage for populations
%
X = zeros( TFinal, 1);  % Initializing storage for population X
tVec = X;               % Initializing storage for time

% Storing initial values
X(1) = x0;
tVec(1) = 0;

%
% For-loop that iteratively solves the discrete dynamical system
%
for n=1:TFinal
   
    % Obtain Next Population
    X(n+1) = X(n) + k*X(n)*( 1 - X(n)/C ); 
    
    % Store Next Time Value
    tVec(n+1) = tVec(n) + 1;
    
end


