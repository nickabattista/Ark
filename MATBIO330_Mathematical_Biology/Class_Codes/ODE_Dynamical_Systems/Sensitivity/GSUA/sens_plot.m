function []=sens_plot(plot_type,p1,p2,p3,p4,p5,p6)
% sens_plot(plot_type,p1,p2,p3,p4,p5,p6) Plots for uncertainty and sensitivity analysis p1,p2,...
% Parametres for every type of plot plot_type
%               sens_plot('FractionalSensitivityArea',Par,SensMethod,S,t)
%               sens_plot('FractionalSensitivityPlots',Par,SensMethod,S,t)
%               sens_plot('TotalSensitivityArea',Par,SensMethod,S,t)
%               sens_plot('TotalSensitivityPlots',Par,SensMethod,S,t)
%               sens_plot('UncertaintyAnalysis',Par,Y,t,y_nom)
%               sens_plot('ScatterParameter',Par,SensMethod,M) - Scatter plot of pair of parameters
%               sens_plot('ScatterOutput',Par,t,SensMethod,Y,M,tref)
%               sens_plot('Bar',Par,SensMethod,ST,t,tref) sens_plot('Pie',Par,SensMethod,ST,t,tref)
% Par           Cell (Npx3) with information about Np parameters
%                   {'Parameter_name', 'Uncertainty_mode', uncertainty_value} Uncertainty_mode = 'range',
%                   uncertainty_value = [min, max] Uncertainty_mode = 'std', uncertainty_value =
%                   [nominal, standard_deviation] Uncertainty_mode = 'percent', uncertainty_value =
%                   [nominal, percent(0,100)]
% SensMethod    Variance-based sensitivity method: 'Saltelli', 'Sobol', 'Jansen', 'Brute-force'
% S             Sensitivity vector (fractional or total): Vi/V (NpxNd)
% Y             Matrix (NdxN) with time responses in columns for every set of parameters
% M             Sample matrix of parameters (one parameter by column and one sample by row):
%                   dim(M)=(NxNp) for brute-force method and dim(M)=(2NxNp) for other methods(M = [A B])
% t             Time vector (Nd)
% tref          Time instant of interest to be plotted
% y_nom         Experimental time response (1xNd) for compute of SJ and STJ. If it is not specified,
%               y_nom is the time response with nominal parameters
%
% (c) Carlos Mario Vélez S. 2015
% Universidad EAFIT, http://www.eafit.edu.co
% Medellín, Antioquia, Colombia
% E-mail: cmvelez@eafit.edu.co
% https://plus.google.com/+CarlosMVelezS/about

Par = p1;
Np = size(Par,1);
switch plot_type
    case 'UncertaintyAnalysis'
        Y = p2; t = p3;
        switch nargin
            case 4
                plot(t,Y,'b');
            case 5
                y_nom = p4;
                h = plot(t,Y,'b',t,y_nom,'r');
                set(h(size(Y,1)+1),'linewidth',2);
            otherwise
                disp('Give the right number of 4 or 5 function inputs')
                disp('sens_plot(''UncertaintyAnalysis'',Par,Y,t)')
                disp('sens_plot(''UncertaintyAnalysis'',Par,Y,t,y_nom)')
                return
        end
        xlabel ('Time')
        title({['Uncertainty analysis (Montecarlo simulation)' ' with N = ' num2str(size(Y,1))]; ' '},'Color','r');
        
    case 'FractionalSensitivityArea'
        if nargin ~=5
            disp('Give the right number of 5 function inputs')
            disp('sens_plot(''FractionalSensitivityArea'',Par,SensMethod,S,t)')
            return
        end
        SensMethod = p2; S = p3; t = p4;
        Nd = size(t,2);
        t0 = 2; % The sensitivity in t = 0 is not computed
        while sum(isnan(S(:,t0))) > 0
            t0 = t0+1;
        end
        area(t(t0:Nd),S(:,t0:Nd)')
        xlabel('Time')
        ylabel('Si = Vi/V')
        title({['Time-dependent fractional sensitivity indices using ' SensMethod ' method']; ' '},'Color','r')
        legend(Par{1:Np},'Location','BestOutside')
        colormap('colorcube')
        
    case 'FractionalSensitivityPlots'
        if nargin ~=5
            disp('Give the right number of 5 function inputs')
            disp('sens_plot(''FractionalSensitivityPlots'',Par,SensMethod,S,t)')
            return
        end
        SensMethod = p2; S = p3; t = p4;
        Nd = size(t,2);
        D1 = floor(sqrt(Np)); % Number of rows of subplot
        D2 = D1+ceil((Np-D1^2)/D1); % Number of columns of subplot
        t0 = 2;
        while sum(isnan(S(:,t0)))>0
            t0 = t0+1;
        end
        for i = 1:Np
            subplot(D1,D2,i)
            plot(t(t0:Nd),S(i,t0:Nd)')
            xlabel('Time')
            ylabel('Si = Vi/V')
            title(Par{i,1})
        end
        h = title(axes,{['Time-dependent fractional sensitivity indices using ' SensMethod ' method'];' '},'Color','r');
        set(gca,'visible','off')
        set(h,'visible','on')
        
    case 'TotalSensitivityArea'
        if nargin ~=5
            disp('Give the right number of 5 function inputs')
            disp('sens_plot(''TotalSensitivityArea'',Par,SensMethod,ST,t)')
            return
        end
        SensMethod = p2; ST = p3; t = p4;
        Nd = size(t,2);
        t0 = 2;
        while sum(isnan(ST(:,t0)))>0
            t0 = t0+1;
        end
        for i = t0:Nd
            ST(:,i)=ST(:,i)./sum(ST(:,i));
        end
        area(t(t0:Nd),ST(:,t0:Nd)')
        xlabel('Time')
        ylabel('SiT = ViT/V   SiT=SiT/sum(SiT)')
        title({['Normalized time-dependent total sensitivity indices using ' SensMethod ' method']; ' '},'Color','r')
        legend(Par{1:Np},'Location','BestOutside')
        colormap('colorcube')
        
    case 'TotalSensitivityPlots'
        if nargin ~=5
            disp('Give the right number of 5 function inputs')
            disp('sens_plot(''TotalSensitivityArea'',Par,SensMethod,ST, t)')
            return
        end
        SensMethod = p2; ST = p3; t = p4;
        D1 = floor(sqrt(Np)); % Number of rows of subplot
        D2 = D1+ceil((Np-D1^2)/D1); % Number of columns of subplot
        for i = 1:Np
            subplot(D1,D2,i)
            plot(t,ST(i,:)')
            xlabel('Time')
            ylabel('SiT = ViT/V')
            title(Par{i,1})
        end
        h = title(axes,{['Time-dependent total sensitivity indices using ' SensMethod ' method'];' '},'Color','r');
        set(gca,'visible','off')
        set(h,'visible','on')
        
    case 'ScatterOutput'  % Scatter plots of Y vs Parameters in the time instant tref
        if nargin ~= 7
            disp('Give the right number of 7 function inputs')
            disp('sens_plot(''ScatterOutput'',Par,SensMethod,Y,M,t,tref)')
            return
        end
        SensMethod = p2; Y = p3; M = p4; t = p5; tref = p6(1);
        D1 = floor(sqrt(Np));
        D2 = D1+ceil((Np-D1^2)/D1);
        [~,c] = find(abs(t-tref)<=min(abs(t-tref))); % Find the column that correspons to a time instant
        for i = 1:Np
            subplot(D1,D2,i)
            plot(M(:,i),Y(:,c),'.')
            xlabel(Par{i,1})
            ylabel('Y')
        end
        h = title(axes,{['Scatterplots of Y vs Parameters in time t = ' num2str(tref) ' with ' SensMethod ' method'];' '},'Color','r');
        set(gca,'visible','off')
        set(h,'visible','on')
        
    case 'ScatterParameter'  % Scatter plots of every pair of parameters
        if nargin ~=4
            disp('Give the right number of 4 function inputs')
            disp('sens_plot(''ScatterParameter'',Par,SensMethod,M)')
            return
        end
        SensMethod = p2; M = p3;
        Ncomb = Np*(Np+1)/2; % Number of combinations of parameters
        D1 = floor(sqrt(Ncomb));
        D2 = D1+ceil((Ncomb-D1^2)/D1);
        l = 1;
        for i = 1:Np
            for j = i:Np
                if j==i
                    subplot(D1,D2,l)
                    hist(M(:,j))
                    xlabel(Par{i,1})
                    ylabel('frequency')
                else
                    subplot(D1,D2,l)
                    plot(M(:,i),M(:,j),'b.')
                    xlabel(Par{i,1})
                    ylabel(Par{j,1})
                end
                l = l+1;
            end
        end
        h = title(axes,{['Scatterplot of pair of parameters with ' SensMethod ' method'];' '},'Color','r');
        set(gca,'visible','off')
        set(h,'visible','on')
        
    case 'Pie'
        switch nargin
            case 6 % Time-dependent pie charts of sensitivity indices
                SensMethod = p2; S = p3; t = p4; tref = p5;
                Nplot = size(tref,2);
                D1 = floor(sqrt(Nplot));
                D2 = D1+ceil((Nplot-D1^2)/D1);
                for j = 1:Nplot
                    % Find the column that correspons to a time instant
                    [~,c] = find( abs(t-tref(j)) <= min(abs(t-tref(j))) );
                    labels = cell(Np);
                    for i = 1:Np
                        if S(i,c)<=0
                            S(i,c) = 1e-6;
                        end
                        if S(i,c)*100/max([sum(S(:,c)),1]) <= 2
                            labels{i}='';
                        else
                            labels{i}=[Par{i} ' (' num2str(S(i,c)*100/max([sum(S(:,c)),1]),3) '%)'];
                        end
                    end
                    subplot(D1,D2,j)
                    labels2 = {'Beta','Alpha','I0'};
                    h1 = pie(S(:,c),labels2);
                    hText = findobj(h1,'Type','text'); % text object handles
                    set(hText,'FontSize',7);
                    title({['Si ( t=' num2str(t(c)) ' )'];' '},'color','red')
                end % Pie chat
                legend(Par{1:Np},'Location','BestOutside')
                h = title(axes,{['Time-dependent pie charts of sensitivity indices using ' SensMethod ' method'];''},'Color','r');
                set(gca,'visible','off')
                set(h,'visible','on')
            case 4 % The pie chart of sensitivity indices for MSE function
                SensMethod = p2; S = p3;
                labels = cell(Np);
                for i = 1:Np
                    if S(i)<=0
                        S(i) = 1e-6;
                    end
                    if S(i)*100/max([sum(S),1]) <= 2
                        labels{i}='';
                    else
                        labels{i}=[Par{i} ' (' num2str(S(i)*100/max([sum(S),1]),3) '%)'];
                    end     
                end
                labels2 = {'Beta','Alpha','I0'};
                labels
                labels={labels{1},labels{2},labels{3}}
                %h1 = pie(S,labels);
                h1=pie(S,labels2);
                hText = findobj(h1,'Type','text'); % text object handles
                set(hText,'FontSize',8);
                title({['Pie chart of sensitivity indices for MSE function by ' SensMethod ' method'];' '},'Color','r')
                legend(Par{1:Np},'Location','BestOutside')
            otherwise
                disp('Give the right number of 4 or 6 function inputs')
                disp('sens_plot(''Pie'',Par,SensMethod,S,t,tref)')
                disp('sens_plot(''Pie'',Par,SensMethod,S)')
                return
        end
        colormap('colorcube')
        
    case 'Bar'
        switch nargin
            case 6 % Time-dependent bar charts of sensitivity indices
                SensMethod = p2; S = p3; t = p4; tref = p5;
                Nplot = size(tref,2);
                D1 = floor(sqrt(Nplot));
                D2 = D1+ceil((Nplot-D1^2)/D1);
                for j = 1:Nplot
                    % Find the column that correspons to a time instant
                    [~,c]  =find(abs(t-tref(j))<=min(abs(t-tref(j))));
                    if sum(S(:,c))<1
                        S1=[S(:,c); 1-sum(S(:,c))];
                    else
                        S1=[S(:,c); 0.00001];
                    end
                    subplot(D1,D2,j)
                    barh(S1)
                    set(gca,'ytick',1:Np,'yticklabel',Par,'FontSize',6)
                    title(['Si ( t=' num2str(t(c)) ' )'])
                    for i = 1:Np
                        if S1(i) > 0.3
                            alignment = 'right';
                            colorl = 'white';
                        else
                            alignment = 'left';
                            colorl = 'blue';
                        end
                        text(S1(i),i,['  ' num2str(S1(i),3)],'HorizontalAlignment',alignment,'color',colorl,'FontSize',6)
                    end
                end
                h = title(axes,{['Time-dependent bar charts of sensitivity indices using ' SensMethod ' method'];' '},'Color','r');
                set(gca,'visible','off')
                set(h,'visible','on')
            case 4 % Bar chart of sensitivity indices for MSE function
                Par = p1; SensMethod = p2; S = p3;
                Np = size(Par,1);
                barh(S);
                set(gca,'ytick',1:Np,'yticklabel',Par,'FontSize',8)
                title({['Bar chart of sensitivity indices for MSE function using ' SensMethod ' method'];' '},'Color','r')
                for i = 1:Np
                    if S(i) > 0.3
                        alignment = 'right';
                        colorl = 'white';
                    else
                        alignment = 'left';
                        colorl = 'blue';
                    end
                    text(S(i),i,['    ' num2str(S(i),3)],'HorizontalAlignment',alignment,'color',colorl,'FontSize',8)
                end
            otherwise
                disp('Give the right number of 4 or 6 function inputs')
                disp('sens_plot(''Bar'',Par,SensMethod,S)')
                disp('sens_plot(''Bar'',Par,SensMethod,S,t,tref)')
                return
        end
        colormap('colorcube')
        
    otherwise
        disp('Select the right plot type')
        return
end
end
