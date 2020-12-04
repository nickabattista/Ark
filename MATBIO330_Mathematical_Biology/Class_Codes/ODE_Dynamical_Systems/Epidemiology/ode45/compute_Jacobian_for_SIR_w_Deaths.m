
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: sets up a Jacobian Matrix and finds eigenvalues
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compute_Jacobian_for_SIR_w_Deaths()


%
% SIR Model Parameters
%
beta = 0.5;                         % rate of disease transmission
gamma = 0.5;                        % rate of recovery
mu = 0.01;                          % natural death rate
muStar = 2*mu;                      % enhanced death rate (natural death rate + death rate from disease)
Lambda = 0.005;                     % birth rate 


%
% Equilibrium Population Values
%
S = ( gamma + muStar ) / beta;
I = Lambda/(gamma+muStar) - mu/beta;
R = (gamma/mu)*I;


%
% Compute elements of Jacobian
%
J11 = -beta*I - mu;             % d(Sdot)/dS
J12 = -beta*S;                  % d(Sdot)/dI
J13 = 0;                        % d(Sdot)/dR
%
J21 = beta*I;                   % d(Idot)/dS
J22 = beta*S - gamma - muStar;  % d(Idot)/dI
J23 = 0;                        % d(Idot)/dR
%
J31 = 0;                        % d(Rdot)/dS
J32 = gamma;                    % d(Rdot)/dI
J33 = -mu;                      % d(Rdot)/dR


%
% Construct Jacobian Matrix
%
J = [J11 J12 J13; J21 J22 J23; J31 J32 J33];


%
% Compute eigenvalues of Jacobian (returns them in a vector array)
%
eigVals = eigs(J);

%
% Print Eigenvalues to Screen
%
fprintf('\n\nTheEIGENVALUES of the JACOBIAN are:');
fprintf('\n\neig1 = %4.4f\n',eigVals(1));
fprintf('\n\neig2 = %4.4f\n',eigVals(2));
fprintf('\n\neig3 = %4.4f\n\n\n',eigVals(3));
%
fprintf('LARGEST eigenvalue is: %4.4f\n\n',max(eigVals));