%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Solves the Logistic Equation using MATLAB's
%           ODE 45 built in differential equation solver, which uses RK-4
%           (4th Order Runge-Kutta Method)
%
%           dP/dt = k*P*(1 - P/C)
%          
%           Parameters: k <- growth rate
%                       C <- carrying capacity
%
% Author: Nick Battista
% Institution: TCNJ
% Created: March 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function go_Go_Logistic_ode45()

%
% Clears any previous plots that are open in MATLAB
clf;

%
% Time Information / Initialization
%
Tstart = 0;            % Simulation starts a tstart (initial value)
Tstop = 150;           % Simulation runs until time = TFinal



%
% Initial Values
%
p0 = 5;                   % Initial Population, P
Initial_Values = [p0];    % Stores initial values in vector



%
% ode45 is matlab's ode solver
%
options=odeset('RelTol',1e-4);
[t,sol] = ode45(@f,[Tstart Tstop],Initial_Values,options);


%
% storing solutions for each variable after solving ODE/ODE System.
% 
P  = sol(:,1);    %gives us P(t)



%
% Plotting solutions
%
plot_Time_Evolutions(t,P)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: RHS vector of the problem: this function evaluates the 
%           RHS of the ODEs and passes it back to the ode45 solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dvdt = f(t,sol)


%
% Components of vectors
%
P = sol(1);         % Susceptible


%
% ODE Parameter Values
%
k = 0.25;           % logistic growth rate
C = 150;            % carrying capacity



%
% ODES (RHS)
%
dPdt = k*P*( 1 - P/C );




%
% Vector to be evaluated
%
dvdt = [dPdt]';



return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the time evolutions (solutions to ODEs)!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Time_Evolutions(t,P)

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
plot(t,P,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Time');
ylabel('Population');
leg = legend('Logistic');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the phase planes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Phase_Planes(S,I,R)

%
% Plot Attributes
%
lw = 4;  % LineWidth (how thick the lines should be)
ms = 25; % MarkerSize (how big the plot points should be)
fs = 18; % FontSize (how big the font should be for labels)

figure(2)
plot(S,I,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Susceptible');
ylabel('Infected');
set(gca,'FontSize',fs);
