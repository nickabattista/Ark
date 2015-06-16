function Eulers()

%Author: Nicholas Battista
%Created: August 12, 2014
%Date of Last Revision: August 23, 2014
%
%This function solves the following ODE:
%dy/dt = f(t,y)
%y(0) = y0
%using Eulers Method. It then does a convergence study for various h values. 
%
%Note: It performs the convergence study for the known ODEs,
%
%dy/dt = y, with y(0) = 1; w/ exact solution is: y(t) = e^t
%and
%dy/dt = 2*pi*cos(2*pi*t), w/ y(0)=1 w/ exact solution y(t) = sin(2*pi*t)+1

print_info();

y0 = 1;   %Initial Condition, y(0)=y0

tS = 0;   %Starting Time for Simulation
tE = 2.0; %End Time for Simulation

%Makes vector of time iterates for convergence study
NVec = give_Me_NumberOfGridPts(); 

%Allocate memory to storage vectors
hVec = zeros(1,length(NVec));
err_h = hVec;
Y_SOL = zeros(1e5,5);
ct = 1;

%For loop for convergence study
for j=1:length(NVec)
    
    %Begin counting time for integration
    tic
    
    N = NVec(j);    %# of time-steps
    h = (tE-tS)/N;  %time-step
    hVec(j) = h;    %stores time-step
    
    time=tS:h:tE;   
    yFOR = zeros(1,length(time));
    yFOR2= yFOR;
    
    %Performs the Euler Method Time-Stepping Scheme
    for i=1:length(time)
    
        if i==1
            yFOR(i) = y0;
            yFOR2(i)= y0;
        else
            yFOR(i) = yFOR(i-1) + h*f(time(i-1),yFOR(i-1),1);
            yFOR2(i)= yFOR2(i-1)+ h*f(time(i-1),yFOR(i-1),2);
        end
        
    end

    %Gives Error at each time-step
    err = compute_Error(time,yFOR,y0,1);
    err2= compute_Error(time,yFOR2,y0,2);

    %Computes Inf-Norm of Error
    err_h(j)  = max( abs(err) );
    err_h2(j) = max( abs(err2) );

    %Stores time for computation
    timeV(j) = toc/2;
    
    %Stores numerical solution / info for plotting soln.
    if ( ( (mod(j,3)==0) && (j<10) ) || (j==11) || (j==14) )
        Y_SOL(1:N+1,ct) = yFOR;
        Y_SOL2(1:N+1,ct)= yFOR2;
        hVecPlot(ct) = h;
        NVecPlot(ct) = N;
        errMat(1:N+1,ct) = err;
        errMat2(1:N+1,ct) = err2;
        ct = ct+1;
    end
    clear h; 
    clear time;
    clear yFOR;

end %Ends looping over different h-values

%Plot Solution vs. Exact
plot_Solution(Y_SOL,Y_SOL2,NVecPlot,hVecPlot,tE,tS,y0,errMat,errMat2)


%Plots Convergence Study
plot_Convergence_Study(hVec,err_h,err_h2);

%Plots time study
plot_Time_Study(NVec,timeV)

fprintf('\nWelp, thats it folks!\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that plots the numerical soln. vs. exact solution for different
% time-step sizes, h, as well as for two different ODEs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Solution(Y_Sol,Y_Sol2,NVec2,hVec2,tE,tS,y0,errMat,errMat2)

fprintf('\nplotting numerical vs. exact solutions...\n\n');

figure(1)
subplot(2,2,1);
strColor = {'k-','-','m-','g-','c-'};
for i=1:length(Y_Sol(1,:))
    h = hVec2(i);
    hVec{i} = strcat('h = ',num2str(h));
    N=  NVec2(i);
    t = tS:h:tE;
    str = strColor{i};
    plot(t,Y_Sol(1:N+1,i),str,'LineWidth',2); hold on;
end
t=tS:0.005:tE;
for i=1:length(t)
    yExact(i) = Exact(t(i),y0,1);
end
plot(t,yExact,'r-','LineWidth',3); hold on;
legend(hVec{1},hVec{2},hVec{3},hVec{4},hVec{5},'Exact Solution','Location','NorthWest');
title('Exact Solns vs Numerical Solns for y(t)=y0*exp(y)');
xlabel('t');
ylabel('y(t)');
%
%
subplot(2,2,2)
for i=1:length(Y_Sol(1,:))
    N=  NVec2(i);
    h = hVec2(i);
    t = tS:h:tE;
    err = abs ( errMat(1:N+1,i) );
    str = strColor{i};
    plot(t,err,str,'LineWidth',2); hold on;
end
legend(hVec{1},hVec{2},hVec{3},hVec{4},hVec{5},'Location','NorthWest');
title('ERROR(t) for various h values for y(t)=y0*exp(y)');
xlabel('t');
ylabel('error(t)');
%
%
%
subplot(2,2,3);
for i=1:length(Y_Sol2(1,:))
    h = hVec2(i);
    N=  NVec2(i);
    t = tS:h:tE;
    str = strColor{i};
    plot(t,Y_Sol2(1:N+1,i),str,'LineWidth',2); hold on;
end
t=tS:0.005:tE;
for i=1:length(t)
    yExact(i) = Exact(t(i),y0,2);
end
plot(t,yExact,'r-','LineWidth',3); hold on;
legend(hVec{1},hVec{2},hVec{3},hVec{4},hVec{5},'Exact Solution','Location','NorthEast');
title('Exact Solns vs Numerical Solns for y(t)=sin(2*pi*t)+1');
xlabel('t');
ylabel('y(t)');
%
%
subplot(2,2,4)
for i=1:length(Y_Sol(1,:))
    N=  NVec2(i);
    h = hVec2(i);
    t = tS:h:tE;
    err = abs ( errMat2(1:N+1,i) );
    str = strColor{i};
    plot(t,err,str,'LineWidth',2); hold on;
end
legend(hVec{1},hVec{2},hVec{3},hVec{4},hVec{5},'Location','NorthEast');
title('ERROR(t) for various h values for y(t)=sin(2*pi*t)+1');
xlabel('t');
ylabel('error(t)');

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that plots the computational time vs. time iterates for [0,T]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Time_Study(Nvec,time)

fprintf('\nplotting computational time study...\n\n');


figure(3)
loglog(Nvec,time,'ro-'); hold on;
xlabel('Number of Time-Steps on [tS,tE]');
ylabel('Log(Time for Each Integration)');
title('Computational Time Study');

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that plots the convergence study, i.e., Error vs. h
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Convergence_Study(hVec,err_h,err_h2)

fprintf('\nplotting convergence studies...\n\n');

figure(2)
subplot(3,4,1)
plot(hVec,err_h,'*'); hold on;
xlabel('h');
ylabel('Inf-Norm Error');
title('Convergence Study: y(t) = exp(t)');

subplot(3,4,2)
semilogx(hVec,err_h,'*'); hold on;
xlabel('Log(h)');
ylabel('Inf-Norm Error');
title('Convergence Study: y(t) = exp(t)');

subplot(3,4,[5,6,9,10])
loglog(hVec,err_h,'*'); hold on;
xlabel('Log(h)');
ylabel('Log(Inf-Norm Error)');
title('Convergence Study: y(t) = exp(t)');

figure(2)
subplot(3,4,3)
plot(hVec,err_h2,'*'); hold on;
xlabel('h');
ylabel('Inf-Norm Error');
title('Convergence Study: y(t) = sin(2*pi*t)+1');

subplot(3,4,4)
semilogx(hVec,err_h2,'*'); hold on;
xlabel('Log(h)');
ylabel('Inf-Norm Error');
title('Convergence Study: y(t) = sin(2*pi*t)+1');

subplot(3,4,[7,8,11,12])
loglog(hVec,err_h2,'*'); hold on;
xlabel('Log(h)');
ylabel('Log(Inf-Norm Error)');
title('Convergence Study: y(t) = sin(2*pi*t)+1');

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that computes Error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function error = compute_Error(time,y,y0,flag)

%Computes Error During Simulation

exact_sol = zeros(1,length(time));
for i=1:length(time)
   exact_sol(i) = Exact(time(i),y0,flag); 
end

error = y - exact_sol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the RHS of the ODE: y' = f(t,y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = f(t,y,flag)

%y' = f(t,y)

if flag == 1
    val = y;
elseif flag==2
    val = 2*pi*cos(2*pi*t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the Exact Sol'n
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = Exact(t,y0,flag)

%Exact sol'n to ODE
if flag == 1
    val = y0*exp(t);
elseif flag == 2
    val = sin(2*pi*t)+1; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function tells you how many time-steps you have for the simulation
% (This is used for the convergence study)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NVec = give_Me_NumberOfGridPts()

N1 = 1:1:9;
N2 = 10:10:90;
N3 = 100:100:900;
N4 = 1000:1000:9*1e3;
N5 = 1e4:1e4:9*1e4;
N6 = 1e5:1e5:9*1e5;

NVec = [N1 N2 N3 N4 N5 N6];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that prints the info about the simulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_info()

fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

fprintf('\nAuthor: Nicholas Battista\n');
fprintf('Created: August 12, 2014\n');
fprintf('Date of Last Revision: August 23, 2014\n\n');

fprintf('This function solves the following ODE:\n\n');
fprintf('dy/dt = f(t,y)\n');
fprintf('y(0) = y0\n\n');
fprintf('using Eulers Method. It then does a convergence study for various h values \n\n');
fprintf('Note: It performs the convergence study for the known ODEs,\n\n');

fprintf('dy/dt = y, with y(0) = 1; so exact solution is: y(t) = e^t\n\n');
fprintf('and\n\n');
fprintf('dy/dt = 2*pi*cos(2*pi*t), w/ y(0)=1 w/ exact solution y(t) = sin(2*pi*t)+1\n\n');

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

