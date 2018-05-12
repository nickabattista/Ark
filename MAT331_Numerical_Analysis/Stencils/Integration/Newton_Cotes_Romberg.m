%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nicholas Battista
% Institution: The College of NJ (TCNJ)
% Email: battistn[at]tcnj[.]edu
% Created: March 25, 2018
%
% This function numerically integrates a function, f(x) between integration
% bounds a and b. This function it integrates is found on line 134.
%
% It computes the Newton-Cotes stencil for a particular
% number of quadrature points, i.e., finds the uniformly spaced quadrature
% pts over [a,b] as well as the quadrature coefficients.
%
% We then compute a composite scheme with just [a,b]->[a,(a+b)/2]U[(a+b)/2,b] 
% and combine the two approximations for Romberg Integration
%
% It will then compare the numerically integrated solution to an exact
% solution if you feed it an exact solution to the integral for testing. It
% then loops over many different Newton-Cotes stencils (i.e., for different
% numbers of quadrature pts and then performs a convergence study to see how
% error decays with higher order stencils.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Newton_Cotes_Romberg()


%Integration Bounds
a = 0;
b = 1;

%Number of Quadrature Pts for Looping
N_S = 1; 
N_E = 18;

%Number of Quadrature Points
err = zeros(1,N_E); err_C = err; err_R = err;
Nvec = N_S:1:N_E;
for j=N_S:N_E-1
    
    %Begin counting time for integration
    tic
    
    %Number of Quad Pts for particular stencil
    N = Nvec(j);
    
    %distance between quad-pts
    dx = (b-a)/(N-1);

    %quad pts
    x = a:dx:b;
    
    %gives vandermond-transpose matrix
    mat = Coeff_Matrix(N,x);

    %gives us RHS to find coeffs
    vec = Monomial_Vector(N,b);

    %gives coefficients
    c = mat\vec;

    %actually does the "integration" on [0,1]
    int = Integrate(N,x,c);

    %computes composite integral [0,1] -> [0,0.5]U[0.5,1]
    x = a:(b-a)/N:b;
    int_composite = Integrate_Composite(N,x,c,2);
    
    %forms Romberg integration approximation
    int_Rom = int + (int_composite - int)/( 1 - (1/2)^N );
    
    %computers error between exact and numerical approximation
    err(N) = abs( (-cos(1)+1) - int );
    err_C(N) = abs( (-cos(1)+1) - int_composite );
    err_R(N) = abs( (-cos(1)+1) - int_Rom );
    
    %stores time for computation
    time(j) = toc;
    
end

% figure attributes
ms = 15; lw = 4; fs = 20;

figure(1)
%subplot(1,4,1:2);
semilogy(Nvec,err(1:end),'o-','MarkerSize',ms,'LineWidth',lw); hold on;
xlabel('Number of Quadrature Pts.');
ylabel('Log(Error)');
title('Convergence Study');
leg = legend('N-pt Stencil');
set(leg,'FontSize',fs-1);
set(gca,'FontSize',fs);
%
pause();
%
semilogy(Nvec(2),err_R(2),'<-','MarkerSize',ms,'LineWidth',lw); hold on;
xlabel('Number of Quadrature Pts.');
ylabel('Log(Error)');
leg = legend('N-pt Stencil','N-pt Stencil w/ Romberg');
title('Convergence Study');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs-1);
%
pause();
%
semilogy(Nvec,err_R(1:end),'<-','MarkerSize',ms,'LineWidth',lw); hold on;
xlabel('Number of Quadrature Pts.');
ylabel('Log(Error)');
leg = legend('N-pt Stencil','N-pt Stencil w/ Romberg');
title('Convergence Study');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs-1);

%figure(2)
%subplot(1,4,3:4);
%semilogy(Nvec,time,'ro-','MarkerSize',ms,'LineWidth',lw); hold on;
%xlabel('Number of Quadrature Pts.');
%ylabel('Log(Time for Each Integration)');
%title('Computational Time Study');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Actually does the Numerical Approximation to the Integral
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int = Integrate(N,x,c)

int = 0;
for i=1:N
   int = int + c(i)*f(x(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Actually does the Numerical Approximation to the Integral using Composite
% Integration with nPartitions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_C = Integrate_Composite(N,x,c,nPartitions)

h = ( x(end) - x(1) )/nPartitions; % Finds length of each sub-interval

int_j = zeros(1,nPartitions);      % Initialize storage fr each sub-integral value

xPartitionPts = x(1):h:x(end);     % Finds partition points in integration domain

for j=1:(nPartitions)
    
    % Make Quadrature Points on Subinterval
    if N==1
        xSub = xPartitionPts(j);
    else
        xSub = xPartitionPts(j):h/(N-1):xPartitionPts(j+1);
    end
    
    % Compute sub-interval integral
    for i=1:N
        int_j(j) = int_j(j) + h*c(i)*f( xSub(i) );
    end
    
end

int_C = sum( int_j );  % Add all sub-interval integral values together to get total composite value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gives us RHS to find coefficients for integration stencil
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vec = Monomial_Vector(N,b)

%This assumes an integration bounds are [0,b]
vec = zeros(N,1);

for i=1:N
   vec(i,1) = b^(i)/i;  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gives us the transpose of vandermonde matrix for finding integration
% stencil coefficients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat = Coeff_Matrix(N,x)

if N==1
    mat = 1;
else

    mat = zeros(N,N);

    for i=1:N
        mat(i,:) = x.^(i-1);  
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to be integrated, f(x)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = f(x) 

val = sin(x);