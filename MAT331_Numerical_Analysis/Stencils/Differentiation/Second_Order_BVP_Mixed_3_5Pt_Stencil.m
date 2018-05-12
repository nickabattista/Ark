%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nick Battista
% Instutition: TCNJ
% Course: MAT 331 (Numerical Analysis)
% Date: 3/4/18
%
% FUNCTION: compute the solution of the following BOUNDARY VALUE PROBLEM:
%           u'' = f
%           u(0) = 1; u(1) = 2;
%
%           NOTE: 1. f is constructed to already know the true solution. 
%                 2. u(x) = (1-x)*cos(x)+2xsin(0.5pi*x)
%                 3. f = 2sin(x) + (x-1)cos(x) + 2pi*cos( 0.5pi*x ) - 0.5pi^2*xsin(0.5pi*x)
%
% Inputs:   N: number of data points to interpolate
%
% Returns:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function maxErr = Second_Order_BVP_Mixed_3_5Pt_Stencil(N)

% Declare boundaries
a=0;
b=1;

% Boundary Values
u0 = 1;
uN = 2;

% Make grid resolution
h = (b-a)/N;

% Construct grid
x = 0:h:1;

% Construct second derivative matrix
D2 = give_Me_Second_Derivative_Matrix(N,h);

% Set up RHS + Boundary Values
RHS = give_Me_RHS(N,x,u0,uN,h);

% Find approximate solution values
uSol = inv(D2)*RHS;

% Tack on EXACT boundary values to ends of approximate solution
uSol = [u0; uSol; uN];

% Plots APPROXIMATE vs. TRUE solution
%plot_Approx_Vs_True_Soln(x,uSol,a,b)

% give MAX ERROR
trueSOL = (1-x).*cos(x)+2*x.*sin(0.5*pi*x);
maxErr = max( abs(trueSOL' - uSol) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plots Approximate vs. True Solution to BVP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Approx_Vs_True_Soln(x,uSol,a,b)

% Plot TRUE solution (on finer grid)
ms = 8; lw = 5;
xFine = a:0.001:b;
plot(xFine,(1-xFine).*cos(xFine)+2*xFine.*sin(0.5*pi*xFine),'r.','MarkerSize',ms,'LineWidth',lw); hold on;

% Plot APPROXIMATE solution
plot(x,uSol,'.','MarkerSize',ms+20,'LineWidth',lw); hold on;
xlabel('x');
ylabel('u(x)');
leg = legend('True Soln','Approx. Soln');
set(gca,'FontSize',18);
set(leg,'FontSize',16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: right hand side of BVP: f(x)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecF = f(x,N)

% Computes RHS of BVP for ALL x
fVec = 2*sin(x) + (x-1).*cos(x) + 2*pi*cos( 0.5*pi*x ) - 0.5*pi^2.*x.*sin(0.5*pi*x);

% Only return the interior of the domain
vecF = fVec(2:N);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Constructs Second Order Derivative Matrix of Size N-1 x N-1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rhs = give_Me_RHS(N,x,u0,uN,h)

% Evaluate RHS of the BVP, f(x).
rhs = f(x,N);

% Tack on the boundary values for left-most and right-most points!
% NOTE: subtraction since coming from LHS of equation
rhs(1) = rhs(1) - u0/h^2;
rhs(2) = rhs(2) - (-1/12)*u0/h^2;
rhs(end-1) = rhs(end-1) - (-1/12)*uN/h^2;
rhs(end) = rhs(end) - uN/h^2;

% Make column vector
rhs = rhs';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Constructs Second Order Derivative Matrix of Size N-1 x N-1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = give_Me_Second_Derivative_Matrix(N,h)

% Initialize A
A = zeros(N-1,N-1);

% Set up tridiagonal matrix
for i=1:N-1
   
    if i==1
        A(i,i) = -2;
        A(i,i+1) = 1;
    elseif i==N-1
        A(i,i) = -2;
        A(i,i-1) = 1;
    elseif i==2
        A(i,i-1) = 4/3; 
        A(i,i) = -5/2;
        A(i,i+1) = 4/3;
        A(i,i+2) = -1/12;
    elseif i==N-2
        A(i,i+1) = 4/3; 
        A(i,i) = -5/2;
        A(i,i-1) = 4/3;
        A(i,i-2) = -1/12;
    else
        A(i,i-2) = -1/12;
        A(i,i-1) = 4/3; 
        A(i,i) = -5/2;
        A(i,i+1) = 4/3;
        A(i,i+2) = -1/12;
    end
        
   
end

% Put factor of 1/h^2 in front
A = 1/h^2*A;
