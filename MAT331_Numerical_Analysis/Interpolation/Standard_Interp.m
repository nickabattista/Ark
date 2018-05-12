%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nick Battista
% Instutition: TCNJ
% Course: MAT 331 (Numerical Analysis)
% Date: 2/5/18
%
% FUNCTION: computes interpolation polynomial using the standard monomial
%           basis
%
% Inputs:   N: number of data points to interpolate
%
% Returns:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Standard_Interp(N)

% initialization
dataX = (0:1/(N-1):1)'; % uniform grid points
dataX = sort(dataX);    % orders vector of x-Values in ascending order
dataY = f(dataX);       % gives corresponding y-data from given function f(x)

% form matrix for finding coefficients (transpose of Vandermonde Matrix)
for i=1:N
    for j=1:N
        A(j,i) = dataX(j)^(i-1);
    end
end

% finds coefficients by inverting the matrix
coeffs = inv(A)*dataY;

% plots TRUE function
lw = 5;
x=0:0.0125:1;
plot( x, f(x), 'b-','LineWidth',lw); hold on; 

pause();

% plots original data points
ms = 10;
plot(dataX,dataY,'k.','MarkerSize',ms+40); hold on;

pause();

% plot INTERPOLATED POLYNOMIAL
x=0:0.0025:1;
for i=1:length(x)
    plot(x(i), p(coeffs,x(i),N),'r.','MarkerSize',ms); hold on;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: function gives function values for input vector, x
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = p(c,x,N)

xPowers = zeros(1,N);
for i=1:N
    xPowers(1,i) = x^(i-1);
end

val = c'*xPowers';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: function gives function values for input vector, x
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = f(x)

val = (x.^2+2).*exp(2*x).*cos(5*x).^2;