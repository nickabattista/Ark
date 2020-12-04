%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Models a Zombie outbreak using Dynamical Systems
%
% Author: Nick Battista
% Created: Jan. 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function go_Go_Zombie_Model(TFinal)

%
% Time Information / Initialization
%TFinal = 100;          % Simulation runs until TFinal
dt = 1e-3;              % Time-Step
TimeVec = 0:dt:TFinal;  % TimeVector = (1,2,3,...,TFinal+1)
H = zeros( TFinal, 1);  % Initializing storage for populations
P = H;
Z = H;

%
% Human birth and death rates
%
b1 = 0.000040;    % Human birth rate (2016)
d1 = 0.000019;    % Natural human death rate (2016)

%
% Human-Zombie Interaction Parameters
%
beta1 = 0.5;  % Human-Zombie Interaction 'Probability' (human -> zombie)
beta2 = 0.05; % Human-Zombie Interaction 'Probability' (human dies, no turning)
beta3 = 0.05; % Human-Zombie Interaction 'Probability' (human kills pre-Zombie)
beta4 = 0.025; % Human-Zombie Interaction 'Probability' (human kills zombie)

%
% pre-Zombie residence time before turning full zombie
%
c1 = 0.1;  % Lag time before pre-Zombie (post-Human) turns Zombie

%
% Zombie death rates
%
d2 = d1*3; % Zombie "natural" death rate


%
% Initial "Populations" (population %)
%
H(1) = 0.99;
P(1) = 0;
Z(1) = 0.01;


for n=1:1:TFinal/dt % TFinal b/c of the way we are array indexing (index 1 = initial time)
   
    % How Human Population Changes
    H(n+1) = H(n) + dt * ( (b1-d1)*H(n) - beta1*H(n)*Z(n) - beta2*H(n)*Z(n) );
    
    % How pre-Zombie Population Changes
    P(n+1) = P(n) + dt * ( - c1*P(n) + beta1*H(n)*Z(n) - beta3*H(n)*P(n) );
    
    % How Zombie Population Changes
    Z(n+1) = Z(n) + dt * ( c1*P(n) - beta4*H(n)*Z(n) - d2*Z(n) );
    
end


%
% PRINT SIMULATION INFO TO SCREEN
%
fprintf('\n\n  --- ZOMBIE MODEL RESULTS --- \n\n');
fprintf('Human Pop. Percentage: %.3f\n\n', 100*H(end) );
fprintf('pre-Zombie Pop. Percentage: %.3f\n\n', 100*P(end) );
fprintf('Zombie Pop. Percentage: %.3f\n\n', 100*Z(end) );
fprintf('Percent Loss (not humans, pre-Z, or Zombies): %.3f\n\n\n', 100*(1-H(end)-P(end)-Z(end) ) );


%
% FIGURE 1: POPULATIONS VS. TIME
%
ms = 30; % MarkerSize for plotting
lw = 4;  % LineWidth for plotting
fs = 18;
%
plot(TimeVec,H,'b.-','MarkerSize',ms,'LineWidth',lw); hold on; 
plot(TimeVec,P,'k.-','MarkerSize',ms,'LineWidth',lw); hold on;
plot(TimeVec,Z,'r.-','MarkerSize',ms,'LineWidth',lw); hold on;
xlabel('Time (hours)');
ylabel('Population %'); 
leg = legend('Humans','pre-Zombies','Zombies');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);


