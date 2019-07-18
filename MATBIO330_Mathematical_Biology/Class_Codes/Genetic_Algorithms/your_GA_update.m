function your_GA_update

clf;                        % clear all figures
N = 5;                      % set the number of simulations to do
fit_val = zeros(1,N);       % create a vector of length N to hold all fitness values
xpos = zeros(1,N);          % create a vector of length N to hold all x positions of the best individual
ypos = zeros(1,N);          % create a vector of length N to hold all y positions of the best individual
savefigs=0;                 % change to 1 to save your figs for each simulation in the same folder as this function

FitnessFunction = @fitness_bird; %this is the name of you function that calcualtes the fitness
numberOfVariables = 2;

%Loop through all N of the simulations.
for i = 1:N
    opts = gaoptimset('PlotFcns',{@gaplotbestf}, 'SelectionFcn',@selectionremainder, 'EliteCount', 5);    %sets your options for the GA
    [xp,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables, opts);                %runs the GA
    Fit_val(i) = Fval;                                                                      %Stores the best fitness at the end of the simulations
    xpos(i) = xp(1);                                                                        %Stores the x values of the best individual    
    ypos(i) = xp(2);                                                                        %Stores the y values of the best individual
    if savefigs==1
        savefig(['Myplot' num2str(i) '.fig']);
    end
end

mean_fit = mean(Fit_val);   %calculate the average best fitness
std_fit = std(Fit_val);     %calculate the standard deviation of the best fitness
min_fit = min(Fit_val);     %determine the most fit individual from all simulations

%print the results in the command window
fprintf('The average best fitness value was : %g\n', mean_fit);
fprintf('The standard deviation of the best fitness value was : %g\n', std_fit);
fprintf('The best fitness value found was : %g\n', min_fit);

%We will use the following information to plot individuals on the fitness graph
M=100;                      %calculate x and y in 100 positions in each direction
x=linspace(0, 10, M);       %hold M values of x between 0 and 10.
y=linspace(0, 10, M);       %hold M values of y between 0 and 10.
z = zeros(M,M);             %create a matrix to hold the fitnesses 

%calculate the fitnesses for each value of x and y.
for i = 1:M
    for j = 1:M
        z(j,i) = (x(i)-5)^2 + (y(j) - 8)^2;
        %z(j,i) = (x(i)-1)^2 * (y(j)-0.5)^2 * (x(i)-0.25)^2 + (y(j)-0.5)^2 + (x(i)-1)^2 * (x(i)-0.25)^2 ;
    end
end


figure(N+1)
hold off
%make a contour plot of the fitness function
contour(x,y,z);
hold on
%plot the x and y positions of the best individuals
plot(xpos,ypos, '*');
xlabel('x');
ylabel('y');
title('Characteristics of best individual in each simulation')


