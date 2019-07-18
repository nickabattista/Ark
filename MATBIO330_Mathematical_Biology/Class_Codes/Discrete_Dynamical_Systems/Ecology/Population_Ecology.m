%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Solves Discrete Dynamical Systems in Population Ecology
%
% Author: Nick Battista
% Institution: TCNJ
% Created: March 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Population_Ecology(TFinal)

%
% Clears any previous plots that are open in MATLAB
clf;

%
% Time Information / Initialization
%
%TFinal = 100;          % Simulation runs until TFinal
TimeVec = 1:1:TFinal;   % TimeVector = (1,2,3,...,TFinal+1)


%
% Initializing storage for populations
%
X = zeros( TFinal, 1);  % Initializing storage for population X
Y = X;                  % Initializing storage for population Y


%
% Initial Values
%
X(1) = 25;
Y(1) = 2;


%
% Parameter Values
%
k = 0.75;    % growth rate
C = 250;    % carrying capacity
b1 = 0.008; % death parameter for prey from predator interactions
b2 = 0.00825; % growth parameter for predator from prey interactions

%
% For-loop that iteratively solves the discrete dynamical system
%
for n=1:TFinal
   
    X(n+1) = X(n) + k*X(n)*( 1 - X(n)/C ) - b1*X(n)*Y(n); 
    Y(n+1) = b2*X(n)*Y(n);
    
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
plot(TimeVec,X(1:end-1),'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
plot(TimeVec,Y(1:end-1),'r.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Time');
ylabel('Population');
leg = legend('Prey','Predator');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);


%
% PLOT 2: Phase Plane Plot
%
figure(2)
plot(X(200:end),Y(200:end),'b.-','LineWidth',lw,'MarkerSize',ms); hold on;
xlabel('Prey Population');
ylabel('Predator Population');
set(gca,'FontSize',fs);
