function Interp()

%Author: Nicholas Battista
%Date of Last Revision: August 14, 2014

%This function interpolates a function, f(x) [on line 228], using the
%Newton Interpolation Scheme. 

%It interpolates over a uniform grid on interval [a,b] for a bunch of
%different studies using an increasing number of interpolation pts. It then
%conducts a convergence study, plotting the error vs. # of interpolation
%pts.

%It will also plot various interpolating polynomials for specific numbers
%of interpolation pts with the data given to interpolate. 


a=-pi/2;  %Starting pt
b= pi/2;  %Ending pt

NS = 1;   %Starting # of Interpolation Pts. - 1 for Conv. Study 
NE = 31;  %Ending # of Interpolation Pts. - 1 for Conv. Study

%Initialization of Storage Matrices / Vectors to use for plotting
C_Mat = zeros(NE,NE);
X_Mat = C_Mat;
err = zeros(NE,1);

for N=NS:NE  % (# of interpolation pts - 1)

    %Uniformly spaced interpolation pts
    h = (b-a)/N;
    x=zeros(1,N);
    x=a:h:b;

    %Associated y-values
    y = f(x);

    %Creates Newton-Interpolation Matrix
    mat =  give_me_Matrix(x);

    %Newton Polynomial Coefficients
    c = mat\y';
    
    %Storing Coefficients / Interpolation Pts
    C_Mat(1:N+1,N) = c;
    X_Mat(1:N+1,N) = x';
    
    %Find Error
    err(N) = compute_Error(a,b,c,x);

end

%Plots Convergence Study
plot_Convergence_Study(NS,NE,err);

%Plots Each Trial Against Actual Poly
plot_Interpolation_Polys(a,b,NS,NE,C_Mat,X_Mat)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to plot the convergence study
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Interpolation_Polys(a,b,NS,NE,C_Mat,X_Mat)

NVec = NS:4:NE; 

len = length(NVec); %# of trials being plotted

%Computes # of subplot columns (tries to make it aesthetically pleasing)
if ceil(len/2) > 4
    Ncols = 3;
else
    Ncols = 2;
end

%Computes # of subplot rows (tries to make it aesthetically pleasing)
Nrows = ceil(len/Ncols);

for j=1:length(NVec); 

    N = NVec(j); %# interpolation pts. - 1     
    c = C_Mat(1:N+1,N);
    x = X_Mat(1:N+1,N);
    h=(b-a)/N;
    
    %Plots Interpolation Pts
    xTest = a:h:b;  
    for i=1:length(xTest)
        yTest(i) =  interp_poly(c,x,xTest(i));
        fTest(i) = f(xTest(i));
    end

    xLen = 0;%(b-a);
    
    %Plots Points along entire interval
    xTest2 = a-xLen/5:h/20:b+xLen/5;
    for i=1:length(xTest2)
        yTest2(i) =  interp_poly(c,x,xTest2(i));
    end
    figure(2)
    subplot(Nrows,Ncols,j);
    plot(xTest,yTest,'.'); hold on;
    plot(xTest,fTest,'r.','MarkerSize',10); hold on;
    plot(xTest2,yTest2,'.','MarkerSize',6); hold on;
    plot(xTest,fTest,'r.','MarkerSize',10); hold on;
    xlabel('x');
    ylabel('Function, f(x) and Poly. fit, p(x)');
    title(sprintf('Poly Fits: N = %d', N))
    legend('p(x)','f(x)');
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to plot the convergence study
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Convergence_Study(NS,NE,err)

subplot(1,2,1);
plot(NS:1:NE,err,'o-');
xlabel('# of pts.');
ylabel('Error');
title('Convergence Study');

subplot(1,2,2);
semilogy(NS:1:NE,err,'o-');
xlabel('# of pts.');
ylabel('log(Error)');
title('Convergence Study');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to evaluate coeff poly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function err = compute_Error(a,b,c,x)

    xTest = a:0.01:b;
    yTest = zeros(1,length(xTest));
    yExact = yTest;
    
    for i=1:length(xTest)
        yTest(i) = interp_poly(c,x,xTest(i));
        yExact(i) = f(xTest(i));
    end
    
    errVec = abs(yTest - yExact);
    err = max(errVec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to evaluate coeff poly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = coeff_poly(i,xV,xPt)

%xPt: pt. we're evaluating the interpolating polynomial at
%xV:  vector of interpolation pts

val = 1;  %initialize
for j=1:i
    val = val * (xPt - xV(j)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to evaluate interpolation polynomial
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function valy = interp_poly(c,xV,xPt)

%xPt: pt. we're evaluating the interpolation polynomial at
%xV: interpolation pts

len = length(c);
valy = c(1);  %initialize
for i=2:len
    valy = valy + c(i)*coeff_poly(i-1,xV,xPt); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to create Matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat =  give_me_Matrix(x)

len = length(x);

mat = ones(len,len);

for i=1:len
    mat(i,2) = x(i) - x(1);
end

for i=3:len    %loops over columns
   for j=1:len %loops over rows 
    
       if (j<i)
        mat(j,i) = 0;
       else
        mat(j,i) = mat(j,i-1)*( x(j) - x(i-1) );
       end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to get data from
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val =  f(x)

val = (x-2).*5.*sin(5*x).*exp(x-1);