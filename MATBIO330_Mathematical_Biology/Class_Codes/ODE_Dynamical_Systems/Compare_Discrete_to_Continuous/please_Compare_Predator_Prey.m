%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Compares Discrete to Continuous Predator-Prey
%
% Author: Nick Battista
% Institution: TCNJ
% Created: March 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function please_Compare_Predator_Prey(b2)


%
% Clears any previous plots that are open in MATLAB
close all;


%
% Time Information / Initialization
%
TFinal = 125;          % Simulation runs until TFinal


%
% Initial Values
%
x0 = 25;
y0 = 2;

%
% Parameter Values
%
k = 0.5;     % growth rate for prey
C = 120;     % carrying capacity for prey
d = 0.1;    % death rate for predator
b1= 0.005;  % prey death rate from interactions w/ predator
%b2= 0.025;  % predator growth rate from interactions w/ prey

%
% Call function to solve Discrete Dynamical System
%
[X_dis,Y_dis,tVec_Discrete] = please_Solve_Discrete_System(TFinal,k,C,d,b1,b2,x0,y0);

%
% Call function to solve Continuous Dynamical System
%
[X_con,Y_con,tVec_Continuous] = please_Solve_Continuous_System(TFinal,k,C,d,b1,b2,x0,y0);

%
% Plot Attributes
%
lw = 3;  % LineWidth (how thick the lines should be)
ms = 30; % MarkerSize (how big the plot points should be)
fs = 18; % FontSize (how big the font should be for labels)

%
% PLOT 1: Prey Populations vs. Time
%
figure(1);
plot(tVec_Discrete,X_dis,'b.','LineWidth',lw,'MarkerSize',ms); hold on;
plot(tVec_Continuous,X_con,'r-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Time');
ylabel('Prey Populations');
leg = legend('Discrete','Continuous');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);

%
% PLOT 2: Predator Populations vs. Time
%
figure(2);
plot(tVec_Discrete,Y_dis,'b.','LineWidth',lw,'MarkerSize',ms); hold on;
plot(tVec_Continuous,Y_con,'r-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Time');
ylabel('Predator Populations');
leg = legend('Discrete','Continuous');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: solves Logistic Differential Equation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X,Y,tVec] = please_Solve_Continuous_System(TFinal,k,C,d,b1,b2,x0,y0)

%
% Initialize Time Information 
%
dt = 0.00125;  % Time-step
t= 0;          % Initial Time

%
% Initialize Initial Population / Time Storage Vector / Counter for While Loop Indexing
%
X(1) = x0;
Y(1) = y0;
tVec(1) = 0;
n = 1;

%
% Solve ODE using Euler's Method
%
while t<TFinal
    
    % Update Time / Counter
    t = t + dt;
    n = n + 1;
    
    %
    % Euler's Method to Solve for Solution
    %
    % Prey
    X(n) = X(n-1) + dt*( k*X(n-1)*( 1 - X(n-1) / C ) - b1*X(n-1)*Y(n-1) ); 
    % Predator
    Y(n) = Y(n-1) + dt*( -d*Y(n-1) + b2*X(n-1)*Y(n-1) );                               
    
    % Update time-Vec
    tVec(n) = tVec(n-1) + dt;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: solves the Discrete Dynamical System for the Logistic Eqn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X,Y,tVec] = please_Solve_Discrete_System(TFinal,k,C,d,b1,b2,x0,y0)

%
% Initializing storage for populations
%
X = zeros( TFinal, 1);  % Initializing storage for prey
Y = X;                  % Initializing storage for predator
tVec = X;               % Initializing storage for time

% Storing initial values
X(1) = x0;
Y(1) = y0;
tVec(1) = 0;

%
% For-loop that iteratively solves the discrete dynamical system
%
for n=1:TFinal
   
    % Obtain Next Population for Prey
    X(n+1) = X(n) + k*X(n)*( 1 - X(n)/C ) - b1*X(n)*Y(n);
    
    % Obtain Next Population for Predator
    Y(n+1) = Y(n) - d*Y(n) + b2*X(n)*Y(n);
    
    % Store Next Time Value
    tVec(n+1) = tVec(n) + 1;
    
end


