%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nick Battista
% Instutition: TCNJ
% Course: MAT 331 (Numerical Analysis)
% Date: 1/20/18
%
% FUNCTION: computes the root of a function using Fixed-Pt Iteration to
%           within a specified error tolerance
%
% Inputs:       p0:   initial guess
%          err_tol:   error tolerance specified    
%             NMAX:   max. number of iterations allowed
%
% Returns:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fixed_Point_Iteration(p0,errTol,NMAX)

% initialization
n = 0;    % counter
err = 1;  % initial error to get started


while ( (n< NMAX) && (err > errTol) )
   
    p = g(p0);              % find next guess at root
    err = abs( p - p0 );    % computes error
    p0 = p;                 % redefine current guess as previous for next loop iteration
    n = n+1;                % add one to counter

end

if n<NMAX
    fprintf('\nThe root is: %d. It took %d iterations\n\n',p0,n);
else
    fprintf('\n Oh no! Fixed Point Iteration Maxed Out...approximate root says: %d...\n\n',p0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: function giving us the fixed-point iteration
%           Examples to solve: x^3 + 4x^2 - 10 = 0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = g(x)

% (b) - doesn't work
%val = sqrt(10/x - 4*x);

% (c)
val = 0.5*sqrt(10-x^3);

% (d)
val = sqrt( 10 / (4+x) );

% (e)
%val = x - (x^3+4*x^2-10) / (3*x^2+8*x);