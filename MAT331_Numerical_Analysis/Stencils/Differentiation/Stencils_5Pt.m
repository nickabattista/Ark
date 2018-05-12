%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nick Battista
% Instutition: TCNJ
% Course: MAT 331 (Numerical Analysis)
% Date: 3/1/18
%
% FUNCTION: computes stencil weights (coefficients) for a 4-point stencil
%           that approximates a function value at xR (the reference point)
%
% Inputs:   
%
% Returns:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Stencils_5Pt()

% 
% SET UP 5-PT STENCIL: f(xR) = c0*f(x0) + c1*f(x1) + c2*f(x2) + c3*f(x3) + c4*f(x4)
% 
% NOTE: using xR = 0.2 (middle point for symmetric stencil)
%

% Stencil Nodes
x0=0;
x1=0.1;
x2=0.2;
x3=0.3;
x4=0.4;

% "Reference" Point
xR=0.10;

% Constuct Vandermonde Transpose Matrix
A1 = [1 1 1 1 1];                     % Row 1
A2 = [x0-xR x1-xR x2-xR x3-xR x4-xR]; % Row 2
A3 = A2.^2;                           % Row 3
A4 = A2.^3;                           % Row 4
A5 = A2.^4;                           % Row 5
A = [A1; A2; A3; A4; A5];  % Pieces the matrix together row by row

% Constructs Right Hand Side, eg, [(x-xR)^0 (x-xR)^1 (x-xR)^2 (x-XR)^3]_{x=0}
b = [1 0 0 0 0]';

% Find coefficients
c = inv(A)*b


% Store coefficints (if you decide to use them for somethng later!)
c0 = c(1);
c1 = c(2);
c2 = c(3);
c4 = c(4);
c5 = c(5);