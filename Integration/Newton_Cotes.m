function Newton_Cotes()

%Author: Nicholas Battista
%Date of Last Revision: August 14, 2014

%This function numerically integrates a function, f(x) between integration
%bounds a and b. This function it integrates is found on line 112.

%It computes the Newton-Cotes stencil for a particular
%number of quadrature points, i.e., finds the uniformly spaced quadrature
%pts over [a,b] as well as the quadrature coefficients.

%It will then compare the numerically integrated solution to an exact
%solution if you feed it an exact solution to the integral for testing. It
%then loops over many different Newton-Cotes stencils (i.e., for different
%numbers of quadrature pts and then performs a convergence study to see how
%error decays with higher order stencils.


%Integration Bounds
a = 0;
b = 1;

%Number of Quadrature Pts for Looping
N_S = 1; 
N_E = 14;

%Number of Quad-pts
err = zeros(1,N_E);
Nvec = N_S:1:N_E;
for j=N_S:N_E
    
    %Begin counting time for integration
    tic
    
    %Number of quadrature pts for particular stencil
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

    %actually does the "integration"
    int = Integrate(N,x,c);
    
    %computers error between exact and numerical approximation
    err(N) = abs( (-cos(1)+1) - int );

    %stores time for computation
    time(j) = toc;
    
end

figure(1)
subplot(1,4,1:2);
semilogy(Nvec,err,'o-'); hold on;
xlabel('Number of Quadrature Pts.');
ylabel('Log(Error)');
title('Convergence Study');

subplot(1,4,3:4);
semilogy(Nvec,time,'ro-'); hold on;
xlabel('Number of Quadrature Pts.');
ylabel('Log(Time for Each Integration)');
title('Computational Time Study');

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

mat = zeros(N,N);

for i=1:N
   mat(i,:) = x.^(i-1);  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to be integrated, f(x)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = f(x) 

val = sin(x);

