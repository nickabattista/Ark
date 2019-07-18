%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script solves the FitzHugh-Nagumo Equations in 1d, which are 
% a simplified version of the more complicated Hodgkin-Huxley Equations. 
%
% Author:  Nick Battista
% Created: 04/21/2019
%
% Equations:
% dv/dt = - v*(v-a)*(v-1) - w + I(t)
% dw/dt = eps*(v-gamma*w)
%
% Variables & Parameters:
% v(t): membrane potential
% w(t): blocking mechanism
% D:      diffusion rate of potential
% a:      threshold potential
% gamma:  resetting rate
% eps:    strength of blocking
% I(t):   initial condition for applied activation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FitzHugh_Nagumo_ODE()

% Parameters in model %
a = 0.15;        % Threshold potential (Note: a=0.3 is traveling wave value, a=0.335 is interesting)
gamma = 1.0;    % Resetting rate (Note: large values give 'funky thick' traveling wave, gamma = 1.0 is desired)
eps = 0.001;    % Blocking strength (Note: eps = 0.001 is desired)
I_mag = 0.05;   % Activation strength

% Temporal  Parameters %
T_final = 5000;      % Sets the final time
Np = 5;              % Set the number of pulses
pulse = T_final/Np;   % determines the length of time between pulses.
NT = 400000;          % Number of total time-steps to be taken
dt = T_final/NT;      % Time-step taken
dp = pulse/50;        % Set the duration of the current pulse
pulse_time = 0;       % pulse time is used to store the time that the next pulse of current will happen

% Initialization %
v = 0;
w = v;
t=0;
tVec = 0:dt:T_final;
Nsteps = length(tVec);
vVec = zeros(Nsteps,1); vVec(1) = v;
wVec = zeros(Nsteps,1); wVec(1) = w;

%
% **** % **** BEGIN SIMULATION! **** % **** %
%
for i=2:Nsteps
    
     % Update the time
    t = t+dt;                        
    
    % Gives applied current activation wave
    [IIapp,pulse_time] = Iapp(pulse_time,I_mag,pulse,dp,t);
    
    % Update potential and blocking mechanism, using Forward Euler
    vN = v + dt * ( - v.*(v-a).*(v-1) - w + IIapp );
    wN = w + dt * ( eps*( v - gamma*w ) );
    
    % Update time-steps
    v = vN;
    w = wN;
    
    % Store time-step values
    vVec(i) = v;
    wVec(i) = w;
    
end

%
% Plots Action Potential (cell voltage) / Blocking Strength vs. Time
%
figure(1)
lw = 5;  % LineWidth
fs = 18; % FontSize
plot(tVec,vVec,'r-','LineWidth',lw); hold on;
plot(tVec,wVec,'b-','LineWidth',lw-2); hold on;
xlabel('Time');
ylabel('Quantity');
leg = legend('Potential','Blocking Strength');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);

%
% Plots Phase Plane: Action Potential vs. Blocking Strength
%
figure(2)
lw = 5;  % LineWidth
fs = 18; % FontSize
plot(wVec,vVec,'k-','LineWidth',lw); hold on;
xlabel('Blocking Strength');
ylabel('Potential');
set(gca,'FontSize',fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: the injection function, Iapp = activation wave for system, and
% returns both the activation as well as updated pulse_time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [app,pulse_time] = Iapp(pulse_time,I_mag,pulse,dp,t)


    %Check to see if there should be a pulse
    if t > (pulse_time)
        
        % Sets pulsing region to current amplitude of I_mag x\in[i1*N,i2*N]
        app = I_mag;  
        
        % Checks if the pulse is over & then resets pulse_time to the next pulse time.
        if t > (pulse_time+dp)
            pulse_time = pulse_time+pulse;
        end
        
    else
        
        % Resets to no activation
        app = 0;
    
    end

