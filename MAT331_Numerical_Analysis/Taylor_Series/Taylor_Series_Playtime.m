%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nick Battista
% Instutition: TCNJ
% Course: MAT 331 (Numerical Analysis)
% Date: 1/10/18
%
% FUNCTION: Finds scaling relation for # of terms in Taylor Series as moves
% away from Maclaurin Centered Point at a=0.
%
% Inputs:        x:   point to evaluate Taylor Series at
%          err_tol: error tolerance specified    
%
% Returns:     num: # of terms necessary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Taylor_Series_Playtime()

% Set an error tolerance   
err_tol = 1e-9;                               

% Decide what x values you want to evaluate 
vals = [1e-2:1e-2:9e-2 1e-1:1e-1:9e-1 1:1:20]; 

% Loop over all possible x values and compute # of terms for specified error
for i=1:length(vals)
    
    % define the x as index i from the vector of possible values
    x = vals(i);      
    
    % use function from earlier to find number of necessary terms in Taylor Series and save it to a new vector
    num_vector(i) = Taylor_Series(x,err_tol);
    
end

%
% plots the x-Value vs # of TS terms
%
figure(1)
plot(vals,num_vector,'*'); hold on;
xlabel('x value put in');
ylabel('number of terms in Taylor Series');

%
% plots the x-Value vs. $ of TS terms but with a semi-log scaling on x-axis
%
figure(2)
semilogx(vals,num_vector,'r*'); hold on;
xlabel('x value put in');
ylabel('number of terms in Taylor Series');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUESTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. What is causing the large jumps in number of terms? (Where it seems discontinuous?
% 2. How do computers actually compute values of sin(100) then?