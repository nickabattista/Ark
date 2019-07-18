%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: solves Diffusion Equation in 2D and compares to Random Walkers
%           in 2D on a lattice
%
%
% Author: Nick Battista
% Institution: TCNJ
% Created: April 8, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Diffusion_2D()


L=10;                         % size of domain [0,L]x[0,L]

ds = 0.05;                    % ds = dx = dy (step-size)

M = 250;                      % # of random walkers

dt = 1e-3;                    % time-step size

x0 = L/2;                     % Initial x-Position of Random Walkers

y0 = L/2;                     % Initial y-Position of Random Walkers

dumpInt = 50;                 % storing dump interval

TFinal = 5;                   % final simulation time

D = ( ds^2 + ds^2 ) / (5*dt); % diffusion coefficient

%
% Solve Diffusion Equation
fprintf('\n\n...Solving 2D Diffusion PDE...\n\n');
[U_store,X,Y,numStored] = please_Solve_Diffusion(L,ds,dt,TFinal,D,dumpInt,M);
fprintf('\n\nFinished solving 2D Diffusion PDE...\n');

pause();

%
% Perform Random Walks
fprintf('\n\n\n...Computing Random Walks...\n\n');
[xMat,yMat,RMS_Sims,RMS_Theory] = please_Perform_Random_Walk(x0,y0,ds,dt,dumpInt,TFinal,M);



%
% Compute Circles of where RMS-Distance is for Simulation
Npts = 100;
for j=2:length(RMS_Sims)
    xC(:,j) =  ( -RMS_Sims(j):2*RMS_Sims(j)/Npts:RMS_Sims(j) )';
    yC_T(:,j) = sqrt( RMS_Sims(j)^2 - xC(:,j).^2 );
    
    xC(:,j) = xC(:,j);
    yC_T(:,j) = yC_T(:,j);
end


%
% Plot Solutions
for j=1:numStored

    uD = U_store(:,:,j);
    
    figure(1)
    
    % Plot Contours for Diffusion and RMS-Simulation Distance Contour
    subplot(1,3,1)
    contourf(X,Y,uD,4,'ShowText','on'); hold on;
    if j==1
        plot(x0,y0,'r.','MarkerSize',20); hold on;
    else
        plot(xC(:,j)+x0,yC_T(:,j)+y0,'r-','LineWidth',6); hold on;
        plot(xC(:,j)+x0,-yC_T(:,j)+y0,'r-','LineWidth',6); hold on;
    end
    axis square;
    
    % Plot Surface in 3D for Diffusion
    subplot(1,3,2)
    surf(X(1,:),Y(:,1),uD); 
    %surf(X,Y,uD,'EdgeColor', 'None', 'facecolor', 'interp'); view(2);
    %axis square;
    
    % Plot Contours for Diffusion and Random Walkers after same amount of time-steps
    subplot(1,3,3)
    contourf(X,Y,uD,4); hold on;
    plot(xMat(:,j),yMat(:,j),'r.','MarkerSize',20); hold on;
    axis square;
    
    % For Class Demonstration
    if j==1
        pause();
    else
        pause(0.5);
    end
    clf;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: perform Random Walks!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xMat,yMat,RMS_SimVec,RMS_Theory] = please_Perform_Random_Walk(x0,y0,ds,dt,dumpInt,TFinal,M)

numSteps = TFinal / dt; % # of steps for Random Walk

ct = 1;              % counter for storing values
xVec = x0*ones(M,1); % initial starting places in x-Position
yVec = y0*ones(M,1); % initial starting places in y-Position
xMat = zeros(M,1);   % initialize storage
yMat = xMat;         % initialize storage
xMat(:,ct) = xVec;   % store initial starting places in x-Position
yMat(:,ct) = yVec;   % store initial starting places in x-Position
RMS_SimVec(ct) = 0;
RMS_Theory(ct) = 0;

for i=1:numSteps
   
    rand_Vec = rand(M);
    
    for j=1:M
       
        % Random Walker Has Choice to Stay in Same Place w/ equal probabiltiy
       if rand_Vec(j) <= 0.20
            xVec(j) = xVec(j) + ds;
       elseif rand_Vec(j) <= 0.4
            xVec(j) = xVec(j) - ds;
       elseif rand_Vec(j) <= 0.6
            yVec(j) = yVec(j) + ds;
       elseif rand_Vec(j) <= 0.8
            yVec(j) = yVec(j) - ds;
       end
        
    end
    
    if mod(i,dumpInt)==0
        ct = ct + 1;
        xMat(:,ct) = xVec;
        yMat(:,ct) = yVec;
        aux = (xVec-x0).^2 + (yVec-y0).^2;
        RMS_SimVec(ct) = sqrt( mean( aux ) );
        RMS_Theory(ct) = sqrt(i)*ds;
        
        fprintf('%d of %d Steps in Random Walk\n',i,numSteps);
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: solve Diffusion PDE in 2D and returns solution and mesh grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U_store,X,Y,ct] = please_Solve_Diffusion(L,ds,dt,TFinal,D,dumpInt,M)

% 
% Give Initial Background Grid for Diffusion
[U,X,Y] = give_Me_Initial_Condition(L,ds,M);

%
% Solve Diffusion Equation
ct = 1;               % counter for storage
U_store(:,:,ct) = U;  % Store initial configuration
t = 0;                % time in simulatio
n = 0;                % number of total time-steps
%
while t<TFinal
   
    n = n+1;
    
    U = Solve_Diffusion(dt,ds,U,D);
    
    if mod(n,dumpInt)==0

        ct = ct+1;

        U_store(:,:,ct) = U;
        
        fprintf('Current Time: %2.3f (of Final Time = %d) for 2D Diffusion\n',t,TFinal);

    end
    
    t = t+dt;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: solves 2D Diffusion Equation w/ Dirichlet Boundary Conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function U = Solve_Diffusion(dt,ds,U,D)

[Ny,Nx] = size(U);

for i=1:Nx

    if i==1

        U(:,1) = 0; %Dirichlet Boundary Condition in x

    elseif i==Nx

        U(:,Nx) = 0; %Dirichlet Boundary Condition in x

    else

        for j=1:Ny

            if j==1

                U(1,:) = 0; %Dirichlet Boundary Condition in y

            elseif j==Ny

                U(Ny,:) = 0; %Dirichlet Boundary Condition in y

            else

                % Solves 2D Diffusion Equation using Finite Differences
                U(j,i) = U(j,i) + dt/ds^2*D*( U(j+1,i) + U(j-1,i) + U(j,i+1) + U(j,i-1) - 4*U(j,i) );

            end

        end

    end
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: give initial condition of background
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,X,Y] = give_Me_Initial_Condition(L,ds,M)

% Give number of grid points
numPts = ceil(L/ds);

% Make Computational Mesh
[X,Y] = meshgrid(0:ds:L,0:ds:L);


%
% Initialize Grid to Zero
U = zeros(numPts+1,numPts+1);

%
% Initial Value of 1 at (x,y)=(L/2,L/2)
U(numPts/2+1,numPts/2+1) = M; 

