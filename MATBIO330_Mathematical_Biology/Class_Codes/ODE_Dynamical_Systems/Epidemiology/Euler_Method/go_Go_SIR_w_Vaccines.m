%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Solves SIR w/ Vaccination (and death) Model
%
% Author: Nick Battista
% Institution: TCNJ
% Created: March 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function go_Go_SIR_w_Vaccines()

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
s0 = 0.99; % Initial Population for Susceptible, S
i0 = 0.1;  % Initial Population for Infected, I
r0 = 0;    % Initial Population for Recovered, R


%
% Parameter Values
%
beta = 0.25;      % rate of disease transmission
gamma = 0.45;     % rate of recovery
nu = 0.25;        % rate of vaccination
mu = 0.0074;      % natural death rate
muStar = 1.5*mu;  % enhanced death rate


%
% Saves Initial Values into Storage Vectors
%
S(1) = s0;
I(1) = i0;
R(1) = r0;
TimeVec(1) = t;

%
% While-loop that iteratively solves the discrete dynamical system
%
while t<TFinal
   
    % Iterate storage counter and time, t
    n = n + 1;
    t = t + dt;
    
    % Define Lambda
    Lambda = mu*S(n) + mu*R(n) + muStar*I(n);
    
    % Solve Susceptible ODE w/ Euler Method
    S(n+1) = S(n) + dt * ( -beta*S(n)*I(n) - mu*S(n) + Lambda -nu*S(n) ); 
    
    % Solve Infected ODE w/ Euler Method
    I(n+1) = I(n) + dt * ( beta*S(n)*I(n) - gamma*I(n) - muStar*I(n) );
    
    % Solve Recovered ODE w/ Euler Method
    R(n+1) = R(n) + dt * ( gamma*I(n) - mu*R(n) + nu*S(n) ); 
    
    % Next time in time storage vector
    TimeVec(n+1) = t;
    
end

%
% Calculate Reproduction Number, R0, for the simulation
%
R0 = beta*Lambda / ((mu+nu)*(gamma+muStar));
fprintf('\n\nReproduction Number, R0 = %4.2f\n\n\n',R0);

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
plot(TimeVec,S,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(TimeVec,I,'r.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(TimeVec,R,'g.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Time');
ylabel('Population');
leg = legend('Susceptible','Infected','Recovered');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);


%
% PLOT 2: Phase Plane Plot
%
figure(2)
plot(S,I,'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Susceptible');
ylabel('Infected');
set(gca,'FontSize',fs);
