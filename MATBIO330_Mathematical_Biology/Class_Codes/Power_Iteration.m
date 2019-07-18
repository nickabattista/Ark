%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes approximations to the largest eigenvalue of a 3x3
%           matrix using the Power Iteration. It then plots those
%           approximations and plots the error betweens successive
%           approximations
%
% Author: Nick Battista
% Institution: TCNJ
% Date Created: February 22, 2019
%
% Inputs: None
%
% Outputs: None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Power_Iteration()

%
% Define Matrix
%
%J = [1 2 3; 4 5 6; 7 8 9];
J = rand(3,3);


%
% MATLAB built in function to compute eigenvalues and then print to screen
%
lambda = eig(J);
fprintf('\nThe eigenvalues of J are: %d, %d, %d\n\n',lambda(1),lambda(2),lambda(3));

%
% Perform Power Iteration
%
v = [1 1 1]';     % initial guess at eigenvector 
err = 1;          % initialize error to get into while-loop
err_Tol = 1e-10;  % error tolerance desired
RQ(1) = 1;        % initial guess at eigenvalue
n = 1;            % counter
while abs(err) > err_Tol
   
    n = n+1;  % counter increases by 1
    
    v = J*v;  % multiply eigenvector guess by Jacobian matrix
    
    RQ(n) = v'*J*v / (v'*v); % compute approximate eigenvalue
    
    err = RQ(n) - RQ(n-1);
    
    err_Vec(n-1) = err;
    
end

%
% Print some info
fprintf('It tooks %d iterations to achieve an error tolerance of %d\n\n',n-1,err_Tol);


%
% plot Dominant Eigenvalue
%
figure(1)
nVec = 1:1:n; % Create vector of iteration numbers
ms = 42;      % MarkerSize for Plotting
lw = 4;       % LineWidth for Plotting
fs = 20;      % FontSize for Plotting
plot(nVec, RQ,'.-','LineWidth',lw,'MarkerSize',ms);
xlabel('Iteration Number');
ylabel('Eigenvalue Estimate');
set(gca,'FontSize',fs);

%
% plot Errors between successive iterations
%
figure(2)
semilogy(nVec(1:end-1), abs(err_Vec),'r.-','LineWidth',lw,'MarkerSize',ms);
xlabel('Iteration Number');
ylabel('Error');
set(gca,'FontSize',fs);

