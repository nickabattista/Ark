%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a non-linear pseudo-spectral ODE solver over all R
%
% Author: Nicholas A. Battista 
% Institution (current): The College of New Jersey (TCNJ)
% Institution (created): Rochester Institute of Technology
% Date Created: August 2009
% Last update: September 2018
%
% Running the code solves the following non-linear ODE:
% 
%      Laplacian(u) + 1/(1+u)^7 = f(r)
%
% with
%
%      f(r) = 8*r^2/(1+r^2)^3 - 6/(1+r^2)^2 + 1/(1+ (1/(1+r^2)) )^7,  [FAST CONV.] 
%      or 
%      f(r) = (r^2-6*r+6)*exp(-r) + 1/(1+(r^2*exp(-r)))^7,  [SLOW CONV.] 
%
% with Dirichelet BCs,
%
%      du/dr(0)=0  &  u(r->inf) = 0 
%
% and exact solution,
%
%      u(r) = 1/(1+r^2)   [FAST CONV.]
%      or
%      u(r) = r^2 e^(-r). [SLOW CONV.]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NonLinear_ODE

Start_Num = 5;  %Number of collocation pts to start simulation
End_Num =  50;  %Number of collocation pts to end simulation

print_info();

%
% Sweeps over different #'s of collocation pts to compare accuracies
for N=Start_Num:End_Num 
    
    %For first iteration, initializes the vectors containing the errors for fast and slow convergence examples
    if N==Start_Num
        NerrorL2_F = zeros(1,N);
        NerrorInf_F = NerrorL2_F;
        time_F = NerrorL2_F;
        
        NerrorL2_S = NerrorL2_F;
        NerrorInf_S = NerrorL2_F;
        time_S = NerrorL2_F;
    end
        
    %Finds solution for particular number of basis functions, N, for FAST convergence example
    fast = 1; %Flag to run the FAST example
    [A un_F NerrorL2_F NerrorInf_F time_F] = find_Solution(N,NerrorL2_F,NerrorInf_F,time_F,fast);
    
    %Finds solution for particular number of basis functions, N, for SLOW convergence example
    fast = 0; %Flag to run the SLOW example
    [A un_S NerrorL2_S NerrorInf_S time_S] = find_Solution(N,NerrorL2_S,NerrorInf_S,time_S,fast);


end %ends for loop at beginning looping over number of grid pts

fprintf('\n -------------------------------------------------------------- \n\n');

plot_collocation_grid(N,A);

plot_solution(N,un_F,un_S)

plot_error_convergence(Start_Num,End_Num,NerrorL2_F,NerrorInf_F,NerrorL2_S,NerrorInf_S);

plot_time_increase(Start_Num,End_Num,time_F,time_S);

plot_coefficients(un_F,un_S)

fprintf('\n\nThat is it! Thanks!\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: finds the approximate solution to the PDF for a particular
%           number of collocation pts, N_A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A un NerrorL2 NerrorInf time] = find_Solution(N_A,NerrorL2,NerrorInf,time,fast)


%collocation points in r and theta directions
A = cheby_collocation_points(N_A);        %%Cheby collocation pts b/w [-1,1]

un = initial_guess(N_A); %Initial guess for spectral coefficients
tol = 1e-8;            %Error Tolerance for Newton's Method
err= 1;                %Error to initialize Newton's Method
n=1;                   %Counter

fprintf('\n -------------------------------------------------------------- \n\n');
fprintf('%d (# of basis functions in x and y)\n\n',N_A);

%Stores Function,Deriv, and 2nd Deriv. Cheby. Poly Values
[TAA TA_P TA_PP] = all_Cheby(A);

if fast == 1
    fprintf('NEWTON METHOD FOR FAST CONVERGENCE EXAMPLE\n');
else
    fprintf('NEWTON METHOD FOR SLOW CONVERGENCE EXAMPLE\n');
end
fprintf('Step | Error\n');

tic
while err > tol
    J = jacobian(N_A,A,TAA,TA_P,TA_PP,un);
    fn = build_rhs(N_A,A,un,fast);
    un1 = un - J\fn;
    err = sqrt((un1-un)'*(un1-un));
    un = un1;

    fprintf('  %d  | %d\n',n,err);
    n=n+1;
end
time(N_A) = toc;

fprintf('Newton Method Converged within tol of %d\n\n',tol);

[NerrorL2 NerrorInf] = expconv(N_A,un,NerrorL2,NerrorInf,fast);

%Nerror = SupNormExpConv(N_A,Nerror,Sinitial,un);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns collocation grid points!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = cheby_collocation_points(N)

x = zeros(1,N+1);
for i=1:N+1
    x(N+2-i) = cos(pi*(i-1)/N);
end
val = x'; %% Need transpose
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the collocation grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_collocation_grid(N_A,A)

fprintf('\nplotting collocation grids...\n');

figure(1)
for k = 1:N_A+1
            plot(A(k),0,'o'); hold on;
end
xlabel('A')
title('Collocation Points')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Builds Jacobian Matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jac = jacobian(N_A,A,TAA,TA_P,TA_PP,un)


jac = zeros(N_A+1,N_A+1);
for a = 1:N_A+1                             %%Runs over Chebyshev Collocation points: A e [-1,1]
             
            for j = 1:N_A+1                 %%Runs over jth Chebyshev Polynomial
                
                TA_pp = TA_PP(j,a);
                TA_p = TA_P(j,a);
                TA = TAA(j,a);
               
                        if a == 1
                              jac(a,j) = TA_p;
                        elseif a == N_A+1
                              jac(a,j) = TA;
                        else
                              %jac(a,j) = 1/4*(1-A(a))^4*TA_pp + 1/2*((1-A(a))^4 /(A(a)+1))*TA_p + 2*interpolateAB(N_A,A(a),un)*TA ;
                              jac(a,j) = 1/4*(1-A(a))^4*TA_pp + 1/2*((1-A(a))^4 /(A(a)+1))*TA_p - 7*(1+interpolateAB(N_A,A(a),un))^(-8)*TA  ;

                        end
                        

            end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Builds Right Hand Side (Boundary Conditions) [ie- Right hand Side of PDE, e.g., Lu = f]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = build_rhs(N_A,A,un,fast)               

rhs = zeros(1,N_A+1);
for a = 1:N_A+1                          %%Runs over Chebyshev Collocation points: A e [0,1]
    
            %UA = exp( - ((1+A(a))/(1-A(a)))^2 );
            %UA_A = 4*(1+A(a))/(A(a)-1)^3*UA;    
            %UA_AA = -8*A(a)*(A(a)^2- 2*A(a) - 7 )/(A(a)-1)^6*UA;
    
    
            if a == 1
                rhs(a) = interpolateAB_p(N_A,A(a),un);
            elseif a == N_A+1
                rhs(a)= interpolateAB(N_A,A(a),un);
            else
                r = (1+A(a))/(1-A(a));
                %rhs(a) = Laplacian(A(a)) + interpolateAB(A(a))^2  - (1/4*(1-A(a))^4*UA_AA + 1/2*((1-A(a))^4/(1+A(a)))*UA_A + UA^2 );
                %rhs(a) = Laplacian(N_A,A(a),un) + interpolateAB(N_A,A(a),un)^2 - (  (r^2-6*r+6)*exp(-r) + (r^2*exp(-r))^2    );
                if fast == 1
                    rhs(a) = Laplacian(N_A,A(a),un) + 1/(1+ interpolateAB(N_A,A(a),un))^7 - ( 8*r^2/(1+r^2)^3 - 6/(1+r^2)^2 + 1/(1+ (1/(1+r^2)) )^7  );    
                else
                    rhs(a) = Laplacian(N_A,A(a),un) + 1/(1+ interpolateAB(N_A,A(a),un))^7 - (  (r^2-6*r+6)*exp(-r) + 1/(1+(r^2*exp(-r)))^7 );
                end
            end
            

end

val = rhs';

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives the value of the approximation series, u(r)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = interpolateAB(N_A,A,un)            

val = 0;

for m =1:N_A+1
    
    TA = T(m-1,A);

    val = val + un(m)*TA;

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives the value of the derivative of the approximation series, u'(r)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = interpolateAB_p(N_A,A,un)            

val = 0;

for m =1:N_A+1
    
    TA_p = Tp(m-1,A);

            val = val + un(m)*TA_p;

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Laplacian in new coordinates, del^2(u)(r)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = Laplacian(N_A,A,un)            

val = 0;

for j = 1:N_A+1
    
            TA_pp = Tpp(j-1,A);
            TA_p = Tp(j-1,A);

            val = val + un(j)*( 1/4*(1-A)^4*TA_pp + 1/2*((1-A)^4 /(A+1))*TA_p  );

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: T''(x) (2nd Derivative of Chebyshev Function of 1st Kind)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = Tpp(j,x)           

if x == 1
    val = (j^4-j^2)/3;
elseif x == -1
    val = (-1)^j*(j^4-j^2)/3;
else
    val = (j*(j+1)*T(j,x)-j*U(j,x))/(x^2-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: T'(x) (1st Derivative of Chebyshev Function of 1st Kind)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = Tp(j,x)           

val = j*U(j-1,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: T(x) (Chebyshev Function of 1st Kind); jth Cheby. polynomial
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = T(j,x)             

if x == 1
    val = 1;
elseif x == -1
    val = (-1)^j;
else
    val = cos(j*acos(x));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: U(x) (Chebyshev Function of 2nd Kind); jth 2nd Cheby. polynomial
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = U(j,x)            

if x == 1
    val = j+1;
elseif x == -1
    val = (-1)^j*(j+1);
elseif j == -1
    val = 0;
else
    val = sin((j+1)*acos(x))/sin(acos(x));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Stores chebyshev polynomials (and 1st & 2nd derivatives) in a
% matrix that is indexed by: 
%                            MAT(k,a): k - which polynomial index, T_k
%                                      a - which collocation pt. index, "X_a"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TA TA_p TA_pp] = all_Cheby(A)

len = length(A);
TA = zeros(len,len);
TA_p = TA;
TA_pp = TA;
for a=1:length(A)
    for k=1:length(A)
        TA(k,a) = T(k-1,A(a));
        TA_p(k,a) = Tp(k-1,A(a));
        TA_pp(k,a) = Tpp(k-1,A(a));
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the approximate solutions vs. the exact solutions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_solution(N,un_F,un_S)

fprintf('\nplotting numerical solns vs. exact soln...\n');

xx = -1:.01:1;

umatrix_F = zeros(1,length(xx));
umatrix_S = umatrix_F;
for i = 1:length(xx)
        umatrix_F(i) = interpolateAB(N,xx(i),un_F);
        umatrix_S(i) = interpolateAB(N,xx(i),un_S);
end

figure(2)
subplot(2,2,1)
plot(xx,umatrix_F,'r')
xlabel('A')
ylabel('U(A)')
title('Numerical Soln: u(r) = 1/(1+r^2) [FAST CONV.]')
axis([-1 1 -0.1 1.1]);

subplot(2,2,2)
ezplot('1/(1+((1+A)/(1-A))^2)',[-1 1])
title('Exact Soln: u(r)=1/(1+r^2) [FAST CONV.]')
axis([-1 1 -0.1 1.1]);

subplot(2,2,3)
plot(xx,umatrix_S,'r')
xlabel('A')
ylabel('U(A)')
title('Numerical Soln: u(r) = r^2 exp(-r) [SLOW CONV.]')
axis([-1 1 -0.1 0.6]);

subplot(2,2,4)
ezplot('((1+A)/(1-A))^2*exp(-((1+A)/(1-A)))',[-1 1])
title('Exact Soln: u(r) = r^2 exp(-r) [SLOW CONV.]')
axis([-1 1 -0.1 0.6]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes the error between approximate solution and exact
%           solution.
%
%           RETURNS: -2 vectors containing L2-error and Inf-error 
%                    -each vector index corresponds to different # of
%                          collocation pts.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NerrorL2, NerrorInf] = expconv(N,un,NerrorL2,NerrorInf,fast)

%Grids to compare solution over
AAA = -1:0.1:1;

oneVector = ones(length(AAA),1);
uAexact = oneVector;

%Computes each separation on variable solution respectively
for i = 1:length(AAA)   
    r = (1+AAA(i))/(1-AAA(i));
    if fast == 1
        uAexact(i) = 1/(1+r^2);   %FAST CONVERGENCE EXAMPLE
    else
        uAexact(i) = r^2*exp(-r); %SLOW CONVERGENCE EXAMPLE
    end
    if i==length(AAA)
       uAexact(i) = 0; 
    end
end

%%%Creates matrix of exact solution @ points [-1,1]x[-1,1] in steps of 0.01
%%%Creates spectral solution @ points [-1,1]x[-1,1] in steps of 0.01
%%%Finds difference between exact and spectral solution
len = length(uAexact);
sol = zeros(1,len);
error = zeros(1,len);
for i = 1:length(AAA)
    sol(i) = interpolateAB(N,AAA(i),un);
    error(i) = abs( sol(i) - uAexact(i) );
end

%Computes L2-Norm Error
absError = sqrt(sum(sum(error)));
fprintf('The L2-Norm Error is: %d\n',absError);

%Computes Inf-Norm Error
maxError = max(max(abs(error)));
fprintf('The Inf-Norm Error is: %d\n',maxError);

%L2-Norm Error Vector
NerrorL2(N) = absError;

%Inf-Norm Error Vector
NerrorInf(N) = maxError;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the Convergence Rates for both cases - "fast" and "slow"
%           convergence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_error_convergence(S,E,NerrorL2_F,NerrorInf_F,NerrorL2_S,NerrorInf_S)

fprintf('\nplotting error convergence...\n');

count = 1:1:E;

figure(3)
%
subplot(2,2,1)
plot(count(S:E),NerrorL2_F(S:E),'r*'); hold on;
plot(count(S:E),NerrorL2_S(S:E),'*');
xlabel('N')
ylabel('L2-Error')
title('Error Convergence: L2-Error vs. N')
legend('fast','slow');
%
subplot(2,2,2)
semilogy(count(S:E),NerrorL2_F(S:E),'r*'); hold on;
semilogy(count(S:E),NerrorL2_S(S:E),'*');
xlabel('N')
ylabel('Log(L2-Error)')
title('Error Convergence: Log(L2-Error) vs. N')
legend('fast','slow');
%
subplot(2,2,3)
plot(count(S:E),NerrorInf_F(S:E),'r*'); hold on;
plot(count(S:E),NerrorInf_S(S:E),'*');
xlabel('N')
ylabel('Inf-Norm Error')
title('Error Convergence: Inf-Norm Error vs. N')
legend('fast','slow');
%
subplot(2,2,4)
semilogy(count(S:E),NerrorInf_F(S:E),'r*'); hold on;
semilogy(count(S:E),NerrorInf_S(S:E),'*');
xlabel('N')
ylabel('Log(Inf-Norm Error)')
title('Error Convergence: Log(Inf-Norm Error) vs. N')
legend('fast','slow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: returns the initial guess for a solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = initial_guess(N_A)

untmp = zeros(1,N_A+1);
for i = 1:N_A+1
    untmp(i) = 0;
end

val = untmp';

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots computational time needed to solve
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_time_increase(S,E,time_F,time_S)

fprintf('\nplotting time complexity...\n');

count = 1:1:E;

figure(4)
subplot(1,2,1);
plot(count(S:E),time_F(S:E),'r*'); hold on;
plot(count(S:E),time_S(S:E),'*'); hold on;
xlabel('N')
ylabel('Time for Each Simulation')
title('Time Complexity vs. N')
legend('fast','slow');
%
subplot(1,2,2);
semilogy(count(S:E),time_F(S:E),'r*'); hold on;
semilogy(count(S:E),time_S(S:E),'*');
xlabel('N')
ylabel('Log(Time for Each Simulation)')
title('Log(Time Complexity) vs. N')
legend('fast','slow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots the spectral coefficients to show exponential decay
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_coefficients(un_F,un_S)

fprintf('\nplotting spectral coefficients...\n');

count = 0:1:length(un_F)-1;


figure(5)
subplot(2,2,1);
plot(count,un_F,'r*'); hold on;
plot(count,un_S,'*'); hold on;
xlabel('n')
ylabel('c_n (coefficient)')
title('Spectral Coefficients: abs(c_n) vs. n')
legend('fast','slow');
axis([0 length(un_F)-1 -0.6 0.6]);

subplot(2,2,2);
plot(count,abs(un_F),'r*'); hold on;
plot(count,abs(un_S),'*'); hold on;
xlabel('n')
ylabel('abs(c_n) (mag. of coefficient)')
title('Spectral Coefficients: abs(c_n) vs. n')
legend('fast','slow');
axis([0 length(un_F)-1 0 0.6]);

subplot(2,2,3:4);
semilogy(count,abs(un_F),'ro'); hold on;
semilogy(count,abs(un_S),'*'); hold on;
xlabel('n')
ylabel('log( abs(c_n) ) (coefficient)')
title('Spectral Coefficients: log(c_n) vs. n')
legend('fast','slow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function print_info()

fprintf('\n\nThis is a non-linear pseudo-spectral ODE solver over all R\n');
fprintf('Author: Nicholas A. Battista \n');
fprintf('Last update: September 2018\n\n');
fprintf('Running the code solves the following non-linear ODE:\n\n');

fprintf('     Laplacian(u) + 1/(1+u)^7 = f(r)\n\n');
fprintf('with\n');
fprintf('     f(r) = 8*r^2/(1+r^2)^3 - 6/(1+r^2)^2 + 1/(1+ (1/(1+r^2)) )^7,  [FAST CONV.] \n\n');
fprintf('     or \n\n');
fprintf('     f(r) = (r^2-6*r+6)*exp(-r) + 1/(1+(r^2*exp(-r)))^7,  [SLOW CONV.] \n\n');
fprintf('with Dirichelet BCs,\n');
fprintf('     du/dr(0)=0, u(r->inf) = 0 \n\n');
fprintf('and exact solution,\n');
fprintf('     u(r) = 1/(1+r^2)  [FAST CONV.]\n');
fprintf('     or \n');
fprintf('     u(r) = r^2 e^(-r)  [SLOW CONV.].\n');
fprintf('\n Note: This simulation will take roughly 5 minutes to complete the convergence study\n');
fprintf('\n -------------------------------------------------------------- \n');
fprintf('\n  -->> BEGIN THE SIMULATION <<--\n');
