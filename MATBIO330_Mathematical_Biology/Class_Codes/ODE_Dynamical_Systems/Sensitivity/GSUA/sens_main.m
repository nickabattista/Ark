clc
fprintf(2,'HELP:\n')
disp('     The Simulink model has to be set correctly (see Configurtion parameters of Simulink examples): ')
disp('     1. The name of parameters are p(1), p(2),...')
disp('     2. Connect an "Out block" to the output')
disp('     3. In "Configuration parameters | Data Import/Export" check these options: ')
disp('        Time (tout), Output (yout), Save simulation output as single object')
disp('     Format of cell of parameters {Parameter1, Uncertainty_mode1, uncertainty_value1; Parameter2, ...}:')
disp('     This tool needs the Statistics Toolbox')
disp(' ')

% Selection of Simulink model and parameters
model_case = input('Select the model name (examples: ''SIR'', ''Pendulum'', ''PID''): ');
switch model_case
    case 'SIR'
        Par = {'beta','percent',[0.15 30];
            'alpha','percent',[0.45 30];
            'I0','percent',[0.1 30]};
        sim_model = 'sens_example_sir_sim';
        y_exp = '';
    case 'Pendulum'
        Par = {'m (kg)','percent',[2.4 20];
            'l (m)','percent',[1.3 20];
            'g (m/s^2)','percent',[9.78,1];
            'f (kg/s)','percent',[0.69,20]};
        sim_model = 'sens_example_pendulum_sim';
        y_exp = '';
    case 'PID'
        Par={'b2','percent',[3.2 20];
            'a1','percent',[2 20];
            'a2','percent',[4 20]};
        sim_model = 'sens_example_pid_sim';
        y_exp = '';
    otherwise
        while exist(model_case,'file')  ~= 4
            fprintf(2,'     The specified Simulink model name does not exit\n')
            model_case = input('Give a correct model name: ');
        end
        Par = input('Give the parameters: ');
        y_exp = input('Give the name of experimental time-response vector ('''' or leave empty for using nominal time response): ');
        sim_model = model_case;
end

% Selection of sensitivity method
SensMethod = input('Select the sensitivity method (''brute-force'', ''Sobol'', ''Jansen'', ''Saltelli''): ');

% Selection of sample method
SampleMethod = input('Select the sample method (''Uniform'', ''LatinHypercube''): ');

% Selection of sample size
N = input('Sample size: ');

% Selection or not of parallel computing
parallel_computing=input('Do use parallel computing? (Yes:1, No: 0): ');

% cluster profile, with the pool size specified by that profile
if parallel_computing == 1 && exist('matlabpool','file')>0 && matlabpool('size')==0
    matlabpool open
end

[S,ST,SJ,STJ,Y,t,M,y_nom] = sens_methods(sim_model,Par,N,SensMethod,SampleMethod,y_exp);

% Plots of relevant information of sensitivity analysis
figure(1)
subplot(2,2,1), sens_plot('UncertaintyAnalysis',Par,Y,t,y_nom)
subplot(2,2,2), sens_plot('FractionalSensitivityArea',Par,SensMethod,S,t)
subplot(2,2,3), sens_plot('TotalSensitivityArea',Par,SensMethod,ST,t)
subplot(2,2,4), sens_plot('Pie',Par,SensMethod,SJ)

figure(2)
sens_plot('ScatterParameter',Par,SensMethod,M)

figure(3)
sens_plot('ScatterOutput',Par,SensMethod,Y,M,t,8)

figure(4)
sens_plot('Pie',Par, SensMethod,ST,t,[1,3,5,7,10])

figure(5)
sens_plot('Bar',Par, SensMethod,ST,t,[1,3,5,7,10])

figure(6)
sens_plot('FractionalSensitivityPlots',Par,SensMethod,S,t)

