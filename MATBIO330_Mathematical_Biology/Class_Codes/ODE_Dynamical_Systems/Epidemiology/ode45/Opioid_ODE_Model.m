%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code solves a model of a basic opioid addiction epidemic
%
%
% Author: Nick Battista
% Date Created: August 10, 2017
% Date Updated: September 23, 2017 (NAB)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Opioid_Basic_Model()

%
% STOCHASTIC? TIME DELAY?
%
global stochastic_flag;
stochastic_flag = 0;
time_delay_flag = 0;

%
% Temporal information
%
Tstart = 0;
Tstop = 100;

%
% initial conditions
%
S_0 = 0.9; 
P_0 = 0.10; 
R_0 = 0.00; 


%
Initial_Values = [S_0 P_0 R_0];


%
% ode45 is matlab's ode solver
%
if time_delay_flag == 0
    options=odeset('RelTol',1e-3);
	[t,sol] = ode45(@f,[Tstart Tstop],Initial_Values,options);
else
    options=odeset('RelTol',1e-3);
    [t,sol] = ode45(@f,[Tstart Tstop],Initial_Values,options); 
end


%
% storing solutions for each variable, theta_k.
% 
S  = sol(:,1);    %gives us S(t)
P =  sol(:,2);    %gives us G(t)
R  = sol(:,3);    %gives us R(t)
A = 1 - S - P - R; %gives us H(t)

%
% Plotting solutions
%
plot_Phase_Planes(S,P,A,R);
plot_Time_Evolutions(t,S,P,A,R)


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: RHS vector of the problem: this function evaluates the 
%           RHS of the ODEs and passes it back to the ode45 solver
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dvdt = f(t,sol)

global stochastic_flag;

%
% ODE Coupling Parameters
%
%
xi = 0.505;
coeff = 0.293;

%
% Dynamical Coupling Parameters
%
alpha = 0.9;        % S->G : people who are prescribed prescription opioids
eps = 0.74;         % G->S : people who use their prescriptions and then go back to susceptible
beta = 0.006;       % S->H : people who get opioids from their relatives/friends/etc to abuse them
mu = 0.00824;       %      : natural death rate
muSTAR = 0.00834 ;  %      : enhanced death rate for opioid abusers
gamma =(1-eps);     % G->H : percent of prescribed opioid class who get addicted to opioids
zeta = 0.75;        % H->R : rate at which Opioid abusers start treatment
delta = 0.09;       % R->S : people who finish their treatment and then go back to susceptible class
nu = coeff*(1-delta);     % R->H : rate at which users in treatment fall back into drug use
sigma = (1-coeff)*(1-delta); % R->H : rate at which people in treatment fall back into use themselves.


% NOTE: sigma+delta+mu = 1.0;
% NOTE: eps+gamma = 1.0;

%
% Components of vectors
%
S = sol(1);         % Susceptible Class
P = sol(2);         % Prescribed Opioid Class
R = sol(3);         % People in Treatment
A = 1 - S - P - R;  % Abusing Opioid Class

Lambda = mu*(S+R+P) + muSTAR*A;

%
% Stochastic Piece ("white noise")
%
if stochastic_flag == 1
    white_noise = awgn(x1,1,'measured');
else
    white_noise = 0;
end


%
% ODES (RHS)
%
Lambda = mu*(S+P+R) + muSTAR*A; 
dS = Lambda - (alpha+mu)*S - beta*(1-xi)*S*A - beta*xi*S*P + eps*P + delta*R;
dP = alpha*S - (gamma+eps+mu)*P;
dR = zeta*A - nu*R*A - (delta+mu+sigma)*R;

% OLD
%dS = Lambda + delta*R - alpha*S - beta*S*H + eps*G - mu*S;
%dP = alpha*S - gamma*G - eps*G - mu*G;
%dR = zeta*H - nu*R*H - delta*R - mu*R - sigma*R;


%
% Vector to be evaluated
%
dvdt = [dS dP dR]';


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots phase planes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Phase_Planes(x1,x2,x3,x4)

figure(1)
plot(x1,x2,'r-','LineWidth',3); hold on;
plot(x1,x3,'b-','LineWidth',3); hold on;
plot(x1,x4,'k-','LineWidth',3); hold on;
xlabel('x1');
ylabel('x2,x3,x4');
legend('G vs. S','H vs. S','R vs. S');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots phase planes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Time_Evolutions(t,x1,x2,x3,x4)

lw = 5;
ms = 10;

figure(2)
plot(t,x1,'-','LineWidth',lw,'MarkerSize',ms,'Color',[0.25    0.5    1]); hold on;
plot(t,x2,'-','LineWidth',lw,'MarkerSize',ms,'Color',[0.9 0.35 0.1]); hold on;
plot(t,x3,'r-','LineWidth',lw,'MarkerSize',ms); hold on; %
plot(t,x4,'-','LineWidth',lw,'MarkerSize',ms,'Color',[0 .4 0]); hold on; %,, ,'Color',[0.3 0 0.9]
xlabel('time');
ylabel('populations');
leg=legend('Susceptible', 'Prescribed', 'Opioid Abuse', 'Treatment');
set(gca,'FontSize',20)
set(leg,'FontSize',18)