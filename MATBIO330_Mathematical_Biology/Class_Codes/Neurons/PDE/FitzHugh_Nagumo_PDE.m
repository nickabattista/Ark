%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script solves the FitzHugh-Nagumo Equations in 1d, which are 
% a simplified version of the more complicated Hodgkin-Huxley Equations. 
%
% Author:  Nick Battista
% Created: 09/11/2015
%
% Equations:
% dv/dt = D*Laplacian(v) - v*(v-a)*(v-1) - w + I(t)
% dw/dt = eps*(v-gamma*w)
%
% Variables & Parameters:
% v(x,t): membrane potential
% w(x,t): blocking mechanism
% D:      diffusion rate of potential
% a:      threshold potential
% gamma:  resetting rate
% eps:    strength of blocking
% I(t):   initial condition for applied activation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FitzHugh_Nagumo_PDE()

% Parameters in model %
D = 1.0;        % Diffusion coefficient
a = 0.3;        % Threshold potential (Note: a=0.3 is traveling wave value, a=0.335 is interesting)
gamma = 1.0;    % Resetting rate (Note: large values give 'funky thick' traveling wave, gamma = 1.0 is desired)
eps = 0.001;    % Blocking strength (Note: eps = 0.001 is desired)
I_mag = 0.05;   % Activation strength

% Discretization/Simulation Parameters %
N = 800;        % # of discretized points   %800
L = 2000;       % Length of domain, [0,L] %500
dx = L/N;       % Spatial Step
x = 0:dx:L;     % Computational Domain

% Temporal  Parameters %
T_final = 10000;      % Sets the final time
Np = 10;              % Set the number of pulses
pulse = T_final/Np;   % determines the length of time between pulses.
NT = 800000;          % Number of total time-steps to be taken
dt = T_final/NT;      % Time-step taken
i1 = 0.475;           % fraction of total length where current starts
i2 = 0.525;           % fraction of total length where current ends
dp = pulse/50;        % Set the duration of the current pulse
pulse_time = 0;       % pulse time is used to store the time that the next pulse of current will happen
IIapp=zeros(1,N+1);   % this vector holds the values of the applied current along the length of the neuron
dptime = T_final/100; % This sets the length of time frames that are saved to make a movie.

% Initialization %
v = zeros(1,N+1);
w = v;
t=0;
ptime = 0;       
tVec = 0:dt:T_final;
Nsteps = length(tVec);
vNext = zeros(Nsteps,N+1); vNext(1,:) = v;
wNext = zeros(Nsteps,N+1); wNext(1,:) = w;

%
% **** % **** BEGIN SIMULATION! **** % **** %
%
for i=2:Nsteps
    
     % Update the time
    t = t+dt;                        
    
    % Give Laplacian
    DD_v_p = give_Me_Laplacian(v,dx);  
    
    % Gives applied current activation wave
    [IIapp,pulse_time] = Iapp(pulse_time,i1,i2,I_mag,N,pulse,dp,t,IIapp);
    
    % Update potential and blocking mechanism, using Forward Euler
    vN = v + dt * ( D*DD_v_p - v.*(v-a).*(v-1) - w + IIapp );
    wN = w + dt * ( eps*( v - gamma*w ) );
    
    % Update time-steps
    v = vN;
    w = wN;
    
    % Store time-step values
    vNext(i,:) = v;
    wNext(i,:) = w;
    
    %This is used to determine if the current time step will be a frame in the movie
    if t > ptime
        figure(1)
        plot(x, v,'r-','LineWidth',5);
        axis([0 L -0.2 1.0]);
        set(gca,'Linewidth',7);
        xlabel('Distance (x)');
        ylabel('Electropotenital (v)');
        ptime = ptime+dptime;
        fprintf('Time(s): %d\n',t);
        pause(0.05);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: the injection function, Iapp = activation wave for system, and
% returns both the activation as well as updated pulse_time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [app,pulse_time] = Iapp(pulse_time,i1,i2,I_mag,N,pulse,dp,t,app)


    %Check to see if there should be a pulse
    if t > (pulse_time)
        
        % Sets pulsing region to current amplitude of I_mag x\in[i1*N,i2*N]
        for j=(floor(i1*N):floor(i2*N))
            app(j) = I_mag;  
        end
        
        % Checks if the pulse is over & then resets pulse_time to the next pulse time.
        if t > (pulse_time+dp)
            pulse_time = pulse_time+pulse;
        end
        
    else
        
        % Resets to no activation
        app = zeros(1,N+1);
    
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives Laplacian of the membrane potential, note: assumes
% periodicity and uses the 2nd order central differencing operator.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DD_v = give_Me_Laplacian(v,dx)

Npts = length(v);
DD_v = zeros(1,Npts);

for i=1:Npts
   if i==1
       DD_v(i) = ( v(i+1) - 2*v(i) + v(end) ) / dx^2;
   elseif i == Npts
       DD_v(i) = ( v(1) - 2*v(i) + v(i-1) ) / dx^2;
   else
       DD_v(i) = ( v(i+1) - 2*v(i) +  v(i-1) ) /dx^2;
   end

end

