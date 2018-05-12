%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nick Battista
% Instutition: TCNJ
% Course: MAT 331 (Numerical Analysis)
% Date: 5/1/18
%
% FUNCTION: computes numerical eigenvalues using QR Algorithm
%
% Inputs:   
%
% Returns:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_QR_Algorithm()

clf; % Clear any existing plots

N = 15;

% CREATE RANDOM MATRIX
A = rand(N,N);
[Q,~] = qr(A); % <-- just to get an orthogonal matrix
v = 10*rand(N,1);
D = diag(v);
A = Q*D*Q';

% SYMMETRIC MATRIX EXAMPLE
%A = randi(2,N,N) - 1;
%A = A - tril(A,-1) + triu(A,1)';

% USE MATLAB TO COMPUTE EIGENVALUES
eigTrue = eig(A);
eigTrue = sort(eigTrue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <<- perform QR Algorithm ->> %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
nIter = 1000;
eigApprox = zeros(N,nIter);

% Loop to Perform QR Algorithm
for j=1:nIter
    [Q,R] = qr(A);
    A = R*Q;
    eigApprox(:,j) = diag(A);
end

% Last Approximation:
eigApprox_Last = eigApprox(:,nIter);
[~,idx] = sort( eigApprox_Last );
eigApprox_Last = eigApprox_Last(idx)

% Error Convergence
[err_L2,err_Inf] = give_Me_Error_For_Every_Iteration(nIter,eigApprox,eigTrue);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <<- PLOTTING THE ERROR ->> %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convergence Plot: L2 Error %%%%
lw = 4; ms = 32; fs = 22;
%
subplot(1,2,1)
plot(1:1:nIter, err_L2,'.-','MarkerSize',ms,'LineWidth',lw); hold on;
xlabel('Iteration Number');
ylabel('L2 Error');
set(gca,'FontSize',fs);
%
subplot(1,2,2)
loglog(1:1:nIter, err_L2,'.-','MarkerSize',ms,'LineWidth',lw); hold on;
xlabel('Iteration Number');
ylabel('L2 Error');
set(gca,'FontSize',fs);

%
% Convergence Plot: L-INFINITY Error %%%%
lw = 4; ms = 32; fs = 22;
%
subplot(1,2,1)
plot(1:1:nIter, err_Inf,'.-','MarkerSize',ms,'LineWidth',lw); hold on;
xlabel('Iteration Number');
ylabel('L2 Error');
set(gca,'FontSize',fs);
leg = legend('L2-Error','L-Inf Error'); 
set(leg,'FontSize',fs);
%
subplot(1,2,2)
loglog(1:1:nIter, err_Inf,'.-','MarkerSize',ms,'LineWidth',lw); hold on;
xlabel('Iteration Number');
ylabel('L2 Error');
set(gca,'FontSize',fs);
leg = legend('L2-Error','L-Inf Error'); 
set(leg,'FontSize',fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: compute the L2 Norm of Error for every iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err_L2,err_Inf] = give_Me_Error_For_Every_Iteration(nIter,eigApprox,eigTrue)

% Initialize L2_Error storage vector
err_L2 = zeros( nIter , 1 );
err_Inf = err_L2;

% ReOrder eigApprox

for j=1:nIter
   
    % Get Approximation Eigenvalues at jth step
    eigApprox_j = eigApprox(:,j);
    
    % Sort Eigenvalues into ascending order (so to compare with MATLAB)
    [~,idx] = sort( eigApprox_j );
    eigApprox_j = eigApprox_j(idx);
    %eigApprox_j = eigApprox_j(end:-1:1);
    
    % Get error between EXACT and APPROX at each step, j
    errVec = eigTrue - eigApprox_j; 
   
    % Compute L2-Error between EXACT AND APPROX at each step, j
    err_L2(j) =  sqrt( errVec'*errVec );
    
    % Compute L-Inf Error between EXACT ANA APPROX at each step, j
    err_Inf(j) = max( abs( errVec ) );
   
end
