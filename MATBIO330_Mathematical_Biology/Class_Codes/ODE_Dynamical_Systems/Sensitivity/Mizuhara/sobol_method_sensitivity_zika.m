%Simulation of Brauer Zika SIR system; Morris sensitivity

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Physical Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%av = [.3, 1]; Gao
%f_hv = [.3, .75]; Gao
%f_vh = [.1, .75]; Gao
%k = [1/12,1/2]; Towers
%g = [1/7,1/3]; Towers
%m = [1/20, 1/6]; Towers
%eta = [1/15, 1/4]; Towers
%a = [.001,.1]; Gao

%bv = av*f_hv = [.09,.75];

    
%b = bv*fvh/f_hv*N/Nv = [.0012, 1.875]; N/Nv ranges from .1 to 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computational parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %Timespan
    Tend = 365;
    tspan = [0 Tend];
    
    N = 10; %number of simulations used in estimate
    
    d = 7; %number of params
    
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sobol Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = sobolset(2*d);

sobol_mat = s(1:N,:); %N by 2d matrix of parameters




    %Rescale
    sobol_mat(:,1)=sobol_mat(:,1)*(1.875-.0012)+.0012; %beta
    sobol_mat(:,2)=sobol_mat(:,2)*(.1-.001)+.001; %alpha
    sobol_mat(:,3)=sobol_mat(:,3)*(1/2-1/12)+1/12;%kappa
    sobol_mat(:,4)=sobol_mat(:,4)*(1/3-1/7)+1/7; %gamma
    sobol_mat(:,5)=sobol_mat(:,5)*(1/6-1/20)+1/20;%mu
    sobol_mat(:,6)=sobol_mat(:,6)*(.75-.09)+.09; %beta_v
    sobol_mat(:,7)=sobol_mat(:,7)*(1/4-1/15)+1/15; %eta
    
    sobol_mat(:,8)=sobol_mat(:,8)*(1.875-.0012)+.0012; %beta
    sobol_mat(:,9)=sobol_mat(:,9)*(.1-.001)+.001; %alpha
    sobol_mat(:,10)=sobol_mat(:,10)*(1/2-1/12)+1/12;%kappa
    sobol_mat(:,11)=sobol_mat(:,11)*(1/3-1/7)+1/7; %gamma
    sobol_mat(:,12)=sobol_mat(:,12)*(1/6-1/20)+1/20;%mu
    sobol_mat(:,13)=sobol_mat(:,13)*(.75-.09)+.09; %beta_v
    sobol_mat(:,14)=sobol_mat(:,14)*(1/4-1/15)+1/15; %eta

   


   tic     
out = sobol_R0(sobol_mat,d,N);

toc

S_R0 = out(1,:);
%ci = out(2:3,:);
S_total_R0=out(2,:);
%ci2 = out(5:6,:);

x=1:d;
figure()
bar(x,S_R0)
%errorbar(x,S_R0,ci(1,:),ci(2,:),'o')
title('R_0: First order Sobol indices')
xticklabels({'\beta','\alpha','\kappa','\gamma','\mu','\beta_v','\eta'})
axis([0 8 0 1.5])
savefig('r0_first_sobol.fig')
saveas(gcf,'r0_first_sobol','png')

figure()
bar(x,S_total_R0)
title('R_0: Total Sobol indices')
xticklabels({'\beta','\alpha','\kappa','\gamma','\mu','\beta_v','\eta'})
axis([0 8 0 1.5])
savefig('r0_total_sobol.fig')
saveas(gcf,'r0_total_sobol','png')

tic
out = sobol_crit_vir(sobol_mat,d,N);
toc

S_R0 = out(1,:);
%ci = out(2:3,:);
S_total_R0=out(2,:);
%ci2 = out(5:6,:);

x=1:d;
figure()
bar(x,S_R0)
title('Critical virulence: First order Sobol indices')
xticklabels({'\beta','\alpha','\kappa','\gamma','\mu','\beta_v','\eta'})
axis([0 8 0 1.5])
savefig('cv_first_sobol.fig')
saveas(gcf,'cv_first_sobol','png')

figure()
bar(x,S_total_R0)
title('Crticial virulence: Total Sobol indices')
xticklabels({'\beta','\alpha','\kappa','\gamma','\mu','\beta_v','\eta'})
axis([0 8 0 1.5])
savefig('cv_total_sobol.fig')
saveas(gcf,'cv_total_sobol','png')


%%%%%%%%%%%%%%%%
%% Sobol Calculations
%%%%%%%%%%%%%%%%
%% Input: A,B, d,N
%% Output: out=[S_y ; S_total]
function out= sobol_y(sobol_mat,d,N)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Matrix A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        A = sobol_mat(1:N,1:d);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Matrix B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        B = sobol_mat(1:N,d+1:2*d);

  mat_size = size(A);

        for i = 1:mat_size(1)

            %Vector of parameters
            X = A(i,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Solve ODE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            y = brauer_zika_ode(X,N,Nv,Tend,I_init);
           % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

           %%%%Measure number of infected humans
            y_A(i) = y(end,3); %Number of infected humans
            % R0_A(i) = X(2)/X(4)+X(1)*X(6)*X(7)/(X(5)*X(4)*(X(5)+X(7)));

        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Matrix B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %B = sobol_mat(1:N,d+1:2*d);
        
            for i = 1:mat_size(1)

            %Vector of parameters
            X = B(i,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Solve ODE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                y = brauer_zika_ode(X,N,Nv,Tend,I_init);
               % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

               %%%%Measure number of infected humans
                y_B(i) = y(end,3); %Number of infected humans
                % R0_B(i) = X(2)/X(4)+X(1)*X(6)*X(7)/(X(5)*X(4)*(X(5)+X(7)));

            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:d %%%each j is a new matrix
    
    A2 = A;
    A2(:,j) = B(:,j); %exchange jth column of A with B
            
    for i = 1:mat_size(1)

            %Vector of parameters
            X = A2(i,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Solve ODE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                y = brauer_zika_ode(X,N,Nv,Tend,I_init);
               % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

               %%%%Measure number of infected humans
                y_mix(i,j) = y(end,3); %Number of infected humans
               %  R0_mix(i,j) = X(2)/X(4)+X(1)*X(6)*X(7)/(X(5)*X(4)*(X(5)+X(7)));

    end

end
    %%y_mix(i,j) : i from 1:N (simulation number) and j from 1:d.
   
    var_y = var(y_A);
    %var_R0 = var(R0_A);
        
    %First order indices
    
    for i=1:d
        
        sumy=0;
        %sumR0 =0;
            for j=1:N

                sumy = sumy+ y_B(j)*(y_mix(j,i)-y_A(j));
               % sumR0 = sumR0+ R0_B(j)*(R0_mix(j,i)-R0_A(j));
            end
        S_y(i) = sumy/N/var_y;
       % S_R0(i) = sumR0/N/var_R0;
    end
    
        
    %Total order indices
    
    for i=1:d
        
        sumy=0;
        %sumR0 =0;
        for j=1:N
        
            sumy = sumy+ (y_A(j)-y_mix(j,i)).^2;
          %  sumR0 = sumR0+ (R0_A(j)-R0_mix(j,i)).^2;
        end
        
        S_total_y(i) = sumy/(2*N)/var_y;
       % S_total_R0(i) = sumR0/(2*N)/var_R0;
    end

    out(1,:)=S_y;
    out(2,:)=S_total_y;
end


%%%%%%%%%%%
%% Use Sobol to calculate R_0 sensitivities
function out=sobol_R0(sobol_mat,d,N)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Matrix A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        A = sobol_mat(1:N,1:d);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Matrix B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        B = sobol_mat(1:N,d+1:2*d);

  mat_size = size(A);

        for i = 1:mat_size(1)

            %Vector of parameters
            X = A(i,:);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Solve ODE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          %  y = brauer_zika_ode(X,N,Nv,Tend,I_init);
           % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

           %%%%Measure number of infected humans
           % y_A(i) = y(end,3); %Number of infected humans
             %X = (beta, alpha,kappa,gamma, mu, beta_v,eta)
             R0_A(i) = R0_calc(X,1,10);
             
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Matrix B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %B = sobol_mat(1:N,d+1:2*d);
        
            for i = 1:mat_size(1)

            %Vector of parameters
            X = B(i,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Solve ODE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               % y = brauer_zika_ode(X,N,Nv,Tend,I_init);
               % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

               %%%%Measure number of infected humans
               % y_B(i) = y(end,3); %Number of infected humans
                R0_B(i) = R0_calc(X,1,10);
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:d %%%each j is a new matrix
    
    A2 = A;
    A2(:,j) = B(:,j); %exchange jth column of A with B
            
    for i = 1:mat_size(1)

            %Vector of parameters
            X = A2(i,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Solve ODE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               % y = brauer_zika_ode(X,N,Nv,Tend,I_init);
               % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

               %%%%Measure number of infected humans
               % y_mix(i,j) = y(end,3); %Number of infected humans
                 R0_mix(i,j) = R0_calc(X,1,10);
    end

end
    %%y_mix(i,j) : i from 1:N (simulation number) and j from 1:d.
   
    %var_y = var(y_A);
    var_R0 = var(R0_A);
        
    %First order indices
    
    %%%%%%%%%%%%%%%%%%%
    %% Matrix of all outputs
    %%%%%%%%%%%%%%%%%%%
    
    final_mat = [R0_A' R0_B' R0_mix];


    %S_R0 = first_order_calc(final_mat,var_R0);
    S_R0 = first_order_ver2(final_mat);
    
    S_total_R0 = total_order_ver2(final_mat);
    
    %%Confidence interval
   % ci = bootci(10000,@(X) first_order_ver2(X),final_mat);
    
   % ci = bootci(1000,@(X) mean_R0(X),final_mat);
    
    %ci2 = bootci(1000,@(X) total_order_calc(X,var_R0),final_mat);
    
    figure()
    histogram(reshape(final_mat,1,N*(d+2)),'Normalization','pdf')
    title('Probability distribution of $R_0$ values','interpreter','latex')
    savefig('r0_pdf.fig')
    saveas(gcf,'r0_pdf','png')

    out(1,:)=S_R0;
   % out(2,:) = ci(1,:);%lower bounds
   % out(3,:) = ci(2,:);   %upper bounds
    out(2,:)=S_total_R0;
    %out(5,:)=ci2(1,:);
    %out(6,:)=ci2(2,:);
end

%%%%%%%%%%%%%%%%%%%%%%
%% Use Sobol to calculate critical virulence
%%%%%%%%%%%%%%%%%%%%%%
function out = sobol_crit_vir(sobol_mat,d,N)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Matrix A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        A = sobol_mat(1:N,1:d);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Matrix B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        B = sobol_mat(1:N,d+1:2*d);

  mat_size = size(A);

        for i = 1:mat_size(1)

            %Vector of parameters
            X = A(i,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Solve ODE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          %  y = brauer_zika_ode(X,N,Nv,Tend,I_init);
           % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

           %%%%Measure number of infected humans
           % y_A(i) = y(end,3); %Number of infected humans
             cv_A(i) = crit_vir_calc(X);
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Matrix B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %B = sobol_mat(1:N,d+1:2*d);
        
            for i = 1:mat_size(1)

            %Vector of parameters
            X = B(i,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Solve ODE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               % y = brauer_zika_ode(X,N,Nv,Tend,I_init);
               % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

               %%%%Measure number of infected humans
               % y_B(i) = y(end,3); %Number of infected humans
                 cv_B(i) = crit_vir_calc(X);
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:d %%%each j is a new matrix
    
    A2 = A;
    A2(:,j) = B(:,j); %exchange jth column of A with B
            
    for i = 1:mat_size(1)

            %Vector of parameters
            X = A2(i,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Solve ODE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               % y = brauer_zika_ode(X,N,Nv,Tend,I_init);
               % [t,f] = ode45(@(t,y) sir_ode(~,y,B1,B2,Bh,b1,b2,bh,nv,nh,n1,n2,m1,m2,mv,mh,dh,a)

               %%%%Measure number of infected humans
               % y_mix(i,j) = y(end,3); %Number of infected humans
                  cv_mix(i,j) = crit_vir_calc(X);
    end

end
    %%y_mix(i,j) : i from 1:N (simulation number) and j from 1:d.
   
    %var_y = var(y_A);
    var_cv = var(cv_A);

    %First order indices
    
    %%%%%%%%%%%%%%%%%%%
    %% Matrix of all outputs
    %%%%%%%%%%%%%%%%%%%
    
    final_mat = [cv_A', cv_B', cv_mix];


    S_R0 = first_order_ver2(final_mat);
    
    S_total_R0 = total_order_ver2(final_mat);
    
    %%Confidence interval
   
    %ci = bootci(10000,{@(X) first_order_calc(X,var_cv),final_mat});
    
    %ci2 = bootci(10000,{@(X) total_order_calc(X,var_cv),final_mat});
    
    figure()
    histogram(reshape(final_mat,1,N*(d+2)),'Normalization','pdf')
    title('Probability distribution of critical virulence values','interpreter','latex')
    savefig('cv_pdf.fig')
    saveas(gcf,'cv_pdf','png')

    
    out(1,:)=S_R0;
    %out(2,:) = ci(1,:);%lower bounds
    %out(3,:) = ci(2,:);   %upper bounds
    out(2,:)=S_total_R0;
    %out(5,:)=ci2(1,:);
    %out(6,:)=ci2(2,:);
end


function out = R0_calc(X,N,Nv)
    %X = (beta, alpha,kappa,gamma, mu, beta_v,eta)
    G = [X(2)/X(4), X(2)/X(4), X(1)*N/Nv*X(7)/(X(5)*(X(5)+X(7))), X(1)*N/Nv/X(5);
        0, 0, 0, 0;
        X(6)*Nv/N/X(4),X(6)*Nv/N/X(4), 0, 0;
        0, 0, 0, 0];
    
    Geigs = eigs(G);
    out = max(Geigs);
end

function [out,index] = crit_vir_calc(X)
    %X = (beta, alpha,kappa,gamma, mu, beta_v,eta)
    
    cycles(1) = X(2)/X(4);
    cycles(2) = sqrt(X(1)*X(6)*X(7)/((X(5)+X(7))*X(4)));
    cycles(3) = (X(1)*X(6)*X(7)/((X(5)+X(7))*X(5)*X(4)))^(1/3);
    
    [out,index]=max(cycles);
    
end


    function out = first_order_calc(final_mat,var_R0)
        size_mat = size(final_mat);
        N = size_mat(1);
        d= size_mat(2)-2;
        out = zeros(d,1);
        
            for i=1:d
        
                %sumy=0;
                sumR0 =0;
                %sumR01=0;
                    for j=1:N

                        %sumy = sumy+ y_B(j)*(y_mix(j,i)-y_A(j));
                        sumR0 = sumR0+ final_mat(j,2)*(final_mat(j,i)-final_mat(j,1));
                        %sumR0 = sumR0+ final_mat(j,2)*final_mat(j,i);
                       % sumR01 = final_mat(j,1);
                    end
                    %sumR0/N^2
                    %var_R0
                   
                %S_y(i) = sumy/N/var_y;
                out(i) = sumR0/N/var_R0;
            end
    
    end
    
    function out = first_order_ver2(final_mat)
        size_mat = size(final_mat);
        Ns = size_mat(1);
       
        for j=1:size_mat(2)-2
            out(j)=(Ns*final_mat(:,2)'*final_mat(:,j+2)-(final_mat(:,2)'*ones(Ns,1))^2)/...
                (Ns*final_mat(:,2)'*final_mat(:,2)-(final_mat(:,2)'*ones(Ns,1))^2);
        end
    end
    
    function out = total_order_ver2(final_mat)
            size_mat = size(final_mat);
        Ns = size_mat(1);
       
        for j=1:size_mat(2)-2
            out(j)=1-...
            (Ns*final_mat(:,1)'*final_mat(:,j+2)-(final_mat(:,2)'*ones(Ns,1))^2)/...
                (Ns*final_mat(:,2)'*final_mat(:,2)-(final_mat(:,2)'*ones(Ns,1))^2);
        end
    end
    
    function out = total_order_calc(final_mat,var_R0)
        size_mat = size(final_mat);
        N = size_mat(1);
        d= size_mat(2)-2;
        
        out = zeros(d,1);
        for i=1:d
        
        %sumy=0;
        sumR0 =0;
        for j=1:N
        
          %  sumy = sumy+ (y_A(j)-y_mix(j,i)).^2;
            sumR0 = sumR0+ (final_mat(j,1)-final_mat(j,i)).^2;
        end
        
       % S_total_y(i) = sumy/(2*N)/var_y;
        out(i) = sumR0/(2*N)/var_R0;
        end
    end