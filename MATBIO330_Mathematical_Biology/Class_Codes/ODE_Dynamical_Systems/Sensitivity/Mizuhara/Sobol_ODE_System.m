%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Performs Sobol Sensitivity Analysis For Calculating Stability
%           of a Disease Free Equilibrium for the SIR Model w/ Deaths
%         
% Orig. Author:  Dr. Matthew Mizuhara (TCNJ)
%
% Modifications: Dr. Nick Battista (TCNJ) 
% Date: April 3, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Sobol_ODE_System()

%
% Computational parameters
%
N = 100000;         % # of simulations used in estimate
d = 5;            % # of parameters used for sensitivity calculations
Tend = 365;       % Time Final (Needed for ODE only)
tspan = [0 Tend]; % Range of Time for ODE to be Solved

%
% Setting Up Sobol Matrix
%

% Constructs a new Sobol sequence point set in d-dimensions.
s = sobolset(2*d);

% Creates N by 2d matrix of parameters
sobol_mat = s(1:N,:); %N by 2d matrix of parameters

% Parameter Order to be Used Below (just for us to reference order here):
% X = [beta gamma mu muStar Lambda]

%
% PARAMETER RANGES
%
% beta
betaLow = 0.05;
betaHigh= 0.5;
% gamma
gammaLow = 0.05;
gammaHigh=0.5;
% mu
muLow = 0.005;
muHigh = 2*muLow;
% muStar 
muStarLow = 2*muLow;
muStarHigh = 10*muLow;
% Lambda
LambdaLow = muLow;
LambdaHigh= 3*muHigh;

%
%Rescale Sobol Matrix Automatically Based on Ranges Above
%
sobol_mat(:,1)=sobol_mat(:,1)*(betaHigh-betaLow)+betaLow;        %beta
sobol_mat(:,2)=sobol_mat(:,2)*(gammaHigh-gammaLow)+gammaLow;     %gamma
sobol_mat(:,3)=sobol_mat(:,3)*(muHigh-muLow)+muLow;              %mu
sobol_mat(:,4)=sobol_mat(:,4)*(muStarHigh-muStarLow)+muStarLow;  %muStar
sobol_mat(:,5)=sobol_mat(:,5)*(LambdaHigh-LambdaLow)+LambdaLow;  %Lambda
%
sobol_mat(:,6)=sobol_mat(:,6)*(betaHigh-betaLow)+betaLow;         %beta
sobol_mat(:,7)=sobol_mat(:,7)*(gammaHigh-gammaLow)+gammaLow;      %gamma
sobol_mat(:,8)=sobol_mat(:,8)*(muHigh-muLow)+muLow;               %mu
sobol_mat(:,9)=sobol_mat(:,9)*(muStarHigh-muStarLow)+muStarLow;   %muStar
sobol_mat(:,10)=sobol_mat(:,10)*(LambdaHigh-LambdaLow)+LambdaLow; %Lambda

   
%
% Actually Perform Sobol Sensitivity (and record how long it takes using tic-toc)
%
fprintf('\n\nTime to run Sobol Sensitivity:\n');
tic     
out = sobol_R0(sobol_mat,d,N);
toc
fprintf('\n\n');


% Store First Order Sobol Indices
S_R0 = out(1,:);

% Store Second Order Sobol Indices
S_total_R0=out(2,:);

% Create Vector of #'s As Dummy for Each Parameter for Plotting
x=1:d;

%
% Creates Bar Graph for 1st Order Sobol Indices
%
figure()
bar(x,S_R0)
title('R_0: First order Sobol indices')
xticklabels({'\beta','\gamma','\mu','\mu^*','\Lambda'})
maxS = 1.05*max(S_R0);
minS = min( 1.05*min(S_R0), 0 );
axis([0 6 minS maxS])

%
% Creates Bar Graph for Total Order Sobol Indices
%
figure()
bar(x,S_total_R0)
title('R_0: Total Sobol indices')
xticklabels({'\beta','\gamma','\mu','\mu^*','\Lambda'})
maxS = 1.05*max(S_total_R0);
minS = min( 1.05*min(S_total_R0), 0 );
axis([0 6 minS maxS])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Performs the Sobol Calculations for ODE System
% 
% Input: A,B, d,N
% Output: out=[S_y ; S_total]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out= sobol_y(sobol_mat,d,N)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A = sobol_mat(1:N,1:d);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    B = sobol_mat(1:N,d+1:2*d);

    mat_size = size(A);

        for i = 1:mat_size(1)

            %Vector of parameters
            X = A(i,:);

            % Solves the ODE System
            y = brauer_zika_ode(X,N,Nv,Tend,I_init);
            % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

            %Measure number of infected humans
            y_A(i) = y(end,3); %Number of infected humans
            % R0_A(i) = X(2)/X(4)+X(1)*X(6)*X(7)/(X(5)*X(4)*(X(5)+X(7)));

        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %B = sobol_mat(1:N,d+1:2*d);
        
        for i = 1:mat_size(1)

            %Vector of parameters
            X = B(i,:);

            % Solves the ODE System
            y = brauer_zika_ode(X,N,Nv,Tend,I_init);
            % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

            % Measure number of infected humans
            y_B(i) = y(end,3); %Number of infected humans
            % R0_B(i) = X(2)/X(4)+X(1)*X(6)*X(7)/(X(5)*X(4)*(X(5)+X(7)));

        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cross matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j=1:d %%%each j is a new matrix

        A2 = A;
        A2(:,j) = B(:,j); %exchange jth column of A with B

        for i = 1:mat_size(1)

            %Vector of parameters
            X = A2(i,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Solve ODE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            y = brauer_zika_ode(X,N,Nv,Tend,I_init);
            % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

            % Measure number of infected humans
            y_mix(i,j) = y(end,3); %Number of infected humans
            %  R0_mix(i,j) = X(2)/X(4)+X(1)*X(6)*X(7)/(X(5)*X(4)*(X(5)+X(7)));

        end

    end
    
    %%y_mix(i,j) : i from 1:N (simulation number) and j from 1:d.
   
    var_y = var(y_A);
    %var_R0 = var(R0_A);
        
    %First order indices
    
    for i=1:d
        
        sumy=0;
        %sumR0 =0;
            for j=1:N

                sumy = sumy+ y_B(j)*(y_mix(j,i)-y_A(j));
               % sumR0 = sumR0+ R0_B(j)*(R0_mix(j,i)-R0_A(j));
            end
        S_y(i) = sumy/N/var_y;
       % S_R0(i) = sumR0/N/var_R0;
    end
    
        
    %Total order indices
    
    for i=1:d
        
        sumy=0;
        %sumR0 =0;
        for j=1:N
        
            sumy = sumy+ (y_A(j)-y_mix(j,i)).^2;
          %  sumR0 = sumR0+ (R0_A(j)-R0_mix(j,i)).^2;
        end
        
        S_total_y(i) = sumy/(2*N)/var_y;
       % S_total_R0(i) = sumR0/(2*N)/var_R0;
    end

    out(1,:)=S_y;
    out(2,:)=S_total_y;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Performs the Sobol Calculations to calculate R_0 sensitivities
% 
% Input: A,B, d,N
% Output: out=[S_y ; S_total]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function out=sobol_R0(sobol_mat,d,N)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A = sobol_mat(1:N,1:d);
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    B = sobol_mat(1:N,d+1:2*d);

    mat_size = size(A);

    for i = 1:mat_size(1)

        %Vector of parameters
        X = A(i,:);

        %
        % Compute R_0: note-> X = (beta, alpha,kappa,gamma, mu, beta_v,eta)
        %        
        R0_A(i) = R0_calc(X);

    end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %B = sobol_mat(1:N,d+1:2*d);
        
    for i = 1:mat_size(1)

        %Vector of parameters
        X = B(i,:);

        %
        % Compute R_0: note-> X = (beta, alpha,kappa,gamma, mu, beta_v,eta)
        %
        R0_B(i) = R0_calc(X);
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cross matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j=1:d %%%each j is a new matrix

        A2 = A;
        A2(:,j) = B(:,j); %exchange jth column of A with B

        for i = 1:mat_size(1)

            %Vector of parameters
            X = A2(i,:);

            %
            % Compute R_0: note-> X = (beta, alpha,kappa,gamma, mu, beta_v,eta)
            %
            R0_mix(i,j) = R0_calc(X);
            
        end

    end
    
    %
    %y_mix(i,j) : i from 1:N (simulation number) and j from 1:d.
    %
    %var_y = var(y_A);
    var_R0 = var(R0_A);
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix of all outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    final_mat = [R0_A' R0_B' R0_mix];

    % First Order Sobol Indices
    S_R0 = first_order_ver2(final_mat);
    
    % Total Sobol Indices
    S_total_R0 = total_order_ver2(final_mat);
    
    %%Confidence interval
    % ci = bootci(10000,@(X) first_order_ver2(X),final_mat);
    % ci = bootci(1000,@(X) mean_R0(X),final_mat);
    % ci2 = bootci(1000,@(X) total_order_calc(X,var_R0),final_mat);
    
    figure()
    histogram(reshape(final_mat,1,N*(d+2)),'Normalization','pdf')
    title('Probability distribution of $R_0$ values','interpreter','latex')

    %
    % Sobol 1st Order and Total Indices
    %
    out(1,:)=S_R0;
    out(2,:)=S_total_R0;
    
    %
    % Confidence Intervals
    %
    % out(2,:) = ci(1,:); %lower bounds
    % out(3,:) = ci(2,:);   %upper bounds
    %out(5,:)=ci2(1,:);
    %out(6,:)=ci2(2,:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns the "R0 value" for a particular parameter set
%
%           Note: not exactly R_0 but a proxy for stability by returning
%                 largest eigenvalue. If eigVal > 0, unstable.
% 
% Input: X,N,Nv
% Output: R_0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R_0 = R0_calc(X)
    
    %
    % SIR w/ Deaths Parameter Vector
    % X = [beta gamma mu muStar Lambda]
    beta = X(1);
    gamma = X(2);
    mu = X(3);
    muStar = X(4);
    Lambda = X(5);
    
    %
    % Disease Free Equilibria (DFE) Pop. Values
    %
    S = Lambda/mu;
    I = 0;
    R = 0;
    
    %
    % Construct Jacobian Matrix
    %
    % 1st Row of Jacobian
    J11 = -beta*I-mu;
    J12 = -beta*S;
    J13 = 0;
    % 2nd Row of Jacobian
    J21 = beta*I;
    J22 = beta*S-gamma-muStar;
    J23 = 0;
    % 3rd Row of Jacobian
    J31 = 0;
    J32 = gamma;
    J33 = -mu;
    
    % Fill in Jacobian Entries
    J = [J11 J12 J13; J21 J22 J23; J31 J32 J33];
    
    % Compute Eigenvalues of Jacobian
    J_eigs = eigs(J);

    % Only Take Largest Eigenvalue (not exactly R_0, but proxy for stability, e.g., epidemic or no)
    R_0 = max(J_eigs);
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes First-Order Sobol Indices
% 
% Input: matrix of all values from trials
% Output: out (first order indices for each parameter)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = first_order_ver2(final_mat)
    
    size_mat = size(final_mat);
    Ns = size_mat(1);

    for j=1:size_mat(2)-2
        out(j)=(Ns*final_mat(:,2)'*final_mat(:,j+2)-(final_mat(:,2)'*ones(Ns,1))^2)/...
            (Ns*final_mat(:,2)'*final_mat(:,2)-(final_mat(:,2)'*ones(Ns,1))^2);
    end
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Total-Order Sobol Indices
% 
% Input: matrix of all values from trials
% Output: out (first order indices for each parameter)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function out = total_order_ver2(final_mat)
            
    size_mat = size(final_mat);
    Ns = size_mat(1);

    for j=1:size_mat(2)-2
        out(j)=1-...
        (Ns*final_mat(:,1)'*final_mat(:,j+2)-(final_mat(:,2)'*ones(Ns,1))^2)/...
            (Ns*final_mat(:,2)'*final_mat(:,2)-(final_mat(:,2)'*ones(Ns,1))^2);
    end
    

    