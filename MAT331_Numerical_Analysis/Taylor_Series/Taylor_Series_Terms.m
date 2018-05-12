%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nick Battista
% Instutition: TCNJ
% Course: MAT 331 (Numerical Analysis)
% Date: 1/10/18
%
% FUNCTION: for a provided error tolerance, computes # of terms necessary
% in Taylor Series for a particular x value
%
% Inputs:        x:   point to evaluate Taylor Series at
%          err_tol: error tolerance specified    
%
% Returns:     num: # of terms necessary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function num = Taylor_Series_Terms(x,err_tol)

% Consider actual function as: sin(x)

% Initialization For While Loop
error = 1;   % could be any value > err_tol
n = -1;      % beginning number of terms in Taylor Series. Note: upon
             %           starting the while loop, will go to n=0

             
% will keep looping and adding an additional term to the Taylor Series until our error is less than that of the error tolerance
while error > err_tol
       
    % Add one more Taylor Series term 
    n = n + 1;

    % Compute Taylor Series
    TS = 0;   %initialize Taylor Series
    for i=0:n
        TS = TS + (-1)^i*(x)^(2*i+1) / factorial(2*i+1);
    end
        
    % Compute the error
    error = abs( sin(x) - TS );
       
end

% If loop ends, last n value must be # of terms until error tolerance was satisfied
num = n;

% Prints to the command window inside of MATLAB
fprintf('It took %d terms to achieve an error tolerance of %d\n',num,err_tol);