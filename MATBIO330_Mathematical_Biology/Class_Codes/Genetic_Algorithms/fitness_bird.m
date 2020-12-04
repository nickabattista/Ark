function [ y ] = fitness_bird( x )

% y is the fitness
% x is a vector of length 3, it holds values for brightness, tail length,
% and song quality.

y = (x(1)-5)^2 + (x(2) - 8)^2;

%y = 4*(x(1)-5)^2 * sin( pi*(x(2)-7) )^4 * (x(2)-13)^2 + 2 * sin(x(1)-5)^2 * ( (x(2)-8)^2 * (x(1) - 5)^2 ) * log( 1 + (x(1)-5)^2 + (x(2)-8)^2*(x(2)-13)^2 ) ; 

%y = (x(1)-1)^2 * (x(2)-0.5)^2 * (x(1)-0.25)^2 + (x(2)-0.5)^2 + (x(1)-1)^2 * (x(1)-0.25)^2;

end

