%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Solves Standard Base Case SIR Model (no deaths) using MATLAB's
%           ODE 45 built in differential equation solver, which uses RK-4
%           (4th Order Runge-Kutta Method)
%
%           dS/dt = Lambda - mu*S - beta*S*I
%           dI/dt = beta*S*I - muStar*S*I - gamma*I
%           dR/dt = gamma*I - mu*R
%          
%           Parameters: Lambda <-- total births added to system (Lambda = mu*S + muStar*I + mu*R)
%                       mu     <-- natural death rate
%                       muStar <-- enhanced death rate (natural death rate + death rate from disease)
%                       beta   <-- rate of disease transmission from interactions btwn healthy and sick person
%                       gamma  <-- rate of recovery
%
% Author: Nick Battista
% Institution: TCNJ
% Created: March 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function go_Go_SIR_ode45()

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
s0 = 0.99;                   % Initial Population for Susceptible, S
i0 = 0.01;                    % Initial Population for Infected, I
r0 = 0;                      % Initial Population for Recovered, R
Initial_Values = [s0 i0 r0]; % Stores initial values in vector



%
% ode45 is matlab's ode solver
%
options=odeset('RelTol',1e-4);
[t,sol] = ode45(@f,[Tstart Tstop],Initial_Values,options);


%
% storing solutions for each variable, theta_k.
% 
S  = sol(:,1);    %gives us S(t)
I =  sol(:,2);    %gives us I(t)
R  = sol(:,3);    %gives us R(t)


%
% Plotting solutions
%
plot_Phase_Planes(S,I,R);
plot_Time_Evolutions(t,S,I,R)




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
S = sol(1);         % Susceptible
I = sol(2);         % Infected
R = sol(3);         % Recovered



%
% ODE Parameter Values
%
beta = 0.55;                        % rate of disease transmission
gamma = 0.45;                       % rate of recovery
mu = 0.0074;                        % natural death rate
muStar = 2*mu;                      % enhanced death rate (natural death rate + death rate from disease)
Lambda = mu*S + muStar*I + mu*R;    % birth rate to equal death rate


%
% ODES (RHS)
%
dSdt = Lambda - mu*S - beta*S*I;
dIdt = beta*S*I - muStar*I - gamma*I;
dRdt = gamma*I - mu*R;



%
% Vector to be evaluated
%
dvdt = [dSdt dIdt dRdt]';



return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the time evolutions (solutions to ODEs)!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Time_Evolutions(t,S,I,R)

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
plot(t,S,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(t,I,'r.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(t,R,'g.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Time');
ylabel('Population');
leg = legend('Susceptible','Infected','Recovered');
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
