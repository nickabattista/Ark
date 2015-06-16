function bingo_bango()

% Author: Nicholas A. Battista
% University: UNC-CH
% First Created: January 19, 2015
% Last Revision: May 5, 2015
%
% This code plays the game of Bingo in a Monte Carlo setting,
% averaging N simulated games of Bingo for a specified number of bingo
% boards.
%
% It tries to answer the following questions:
% 1. What is the expected number of draws to win regular bingo if there are
% N bingo boards in play?
% 2. What is the convergence rate of expected number of draws vs. number of
% bingo boards for a certain number of simulated bingo games?
% 3. What about inner square, outer square, and cover all?
%

%print simulation information
print_info();


BoardsVec = [1 5:5:75];       %Set the number of bingo boards for simulation
NgamesVec = [5 10 25 50 100]; %Set the number of games in each series of simulations


for k=1:length(NgamesVec)

    Ngames = NgamesVec(k);
    
    fprintf('\n\n-----------------------------------------------------------------------------------\n');
    fprintf('\n******** Playing %d games of bingo for each number of bingo cards (players) ******** ',Ngames);
    fprintf('\n\n-----------------------------------------------------------------------------------\n\n');

    for i=1:length(BoardsVec)
    
        Nboards = BoardsVec(i);
        
        if Nboards == 1
            fprintf('-> Playing with %d bingo cards\n',Nboards);
        else
            fprintf('-> Playing with %d bingo cards\n',Nboards);
        end
        
        for j=1:Ngames

            B = please_Make_Boards(Nboards); %Makes boards for game!
            balls = randperm(75);            %Draw balls randomly for entire game!
            
 
            %Play Regular Bingo -> Return boards and counter for inner/outer square/cover-all
            [B ct_reg num_winners_reg] = play_Regular_Bingo(B,balls);
            
            [B ct_inner num_winners_in] = play_Inner_Square_Bingo(B,balls,ct_reg);
            
            [B ct_outer num_winners_out] = play_Outer_Square_Bingo(B,balls,ct_inner);
            
            [B ct_cover num_winners_cover] = play_Cover_All_Bingo(B,balls,ct_outer);
            
            
            %Row: different games for specific # of boards
            %Column: Different # of boards
            %k: Different amount of games played in a series of simulations for convergence
            count_reg(j,i,k) = ct_reg;
            numWinnersReg(j,i,k) = num_winners_reg;

            count_inner(j,i,k) = ct_inner;
            numWinnersInner(j,i,k) = num_winners_in;
            
            count_outer(j,i,k) = ct_outer;
            numWinnersOut(j,i,k) = num_winners_out;
            
            count_cover(j,i,k) = ct_cover;
            numWinnersCover(j,i,k) = num_winners_cover;

        end
    end
end


% AvgMat: row - different # of Monte carlo sims. per game / column: different # of bingo cards
AvgMatReg = please_Compute_Averages(count_reg,NgamesVec);
AvgMatIn = please_Compute_Averages(count_inner,NgamesVec);
AvgMatOut = please_Compute_Averages(count_outer,NgamesVec);
AvgMatCover = please_Compute_Averages(count_cover,NgamesVec);

% Plots data of interest
plot_Games_Required_To_Win(NgamesVec,BoardsVec,AvgMatReg,AvgMatIn,AvgMatOut,AvgMatCover);
plot_Convergence_Study(BoardsVec,AvgMatReg,AvgMatIn,AvgMatOut,AvgMatCover);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: constructs Nboards of bingo randomly and stores them in a 3D
% Matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = please_Make_Boards(Nboards)

%Nboards: Number of boards

B = zeros(5,5,Nboards);

%B: 1-15
%I: 16-30
%N: 31-45
%G: 46-60
%O: 61-75

for i=1:Nboards

    %Creates random vector of numbers in appropriate interval
    B1 = randperm(15);
    I1 = randperm(15)+15;
    N1 = randperm(15)+30;
    G1 = randperm(15)+45;
    O1 = randperm(15)+60;

    %Take first 5 random values in appropriate intervals for board columns
    Bvec = B1(1:5)';
    Ivec = I1(1:5)';
    Nvec = N1(1:5)'; 
    Gvec = G1(1:5)';
    Ovec = O1(1:5)';

    %Set FREE SPACE (to zero)
    Nvec(3) = 0;

    %Sets the Bingo Board!
    B(:,1,i) = Bvec;
    B(:,2,i) = Ivec;
    B(:,3,i) = Nvec;
    B(:,4,i) = Gvec;
    B(:,5,i) = Ovec;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Plays "regular" bingo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B ct num_winners] = play_Regular_Bingo(B,balls)

win_flag = 0;         %initialize winning flag to end game
ct = 0;               %counter for number of balls drawn

while ( win_flag == 0 )

    ct = ct+1;        %updates counter when picking another ball
    ball = balls(ct); %draws a ball from already randomly drawn balls
    B = please_Dab_Bingo_Boards(B,ball); %dabs bingo board (sets elements to ZERO)
        
    %Checks for Regular Bingo
    [win_flag,num_winners] = please_Check_Regular_Bingo_For_Win(B); 
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCION: "dabs" out bingo board when certain ball is called. 
% Note: Board starts with all 1's, Change 1 -> 0 when number is called.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = please_Dab_Bingo_Boards(B,ball)

indices = find(B==ball);
B(indices) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Checks to see if any boards have won regular bingo.
% Regular Bingo: Straight Lines (horizontal, vertical, diagonal), inner
% and/or outer four corners.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [win_flag,num_winners] = please_Check_Regular_Bingo_For_Win(B)

Nboards = length(B(1,1,:)); 
win_flag = 0;
num_winners = 0;
win_flag2 = 0;

for i=1:Nboards
    
    %Check ROWS
    if     max(B(1,:,i))==0
        win_flag2 = 1;
    elseif max(B(2,:,i))==0
        win_flag2 = 1;
    elseif max(B(3,:,i))==0
        win_flag2 = 1;
    elseif max(B(4,:,i))==0
        win_flag2 = 1;
    elseif max(B(5,:,i))==0
        win_flag2 = 1;
    end
    
    %Check COLUMNS if board hasn't won yet    
    if win_flag2 == 0
        if     max(B(:,1,i))==0
            win_flag2 = 1;
        elseif max(B(:,2,i))==0
            win_flag2 = 1;
        elseif max(B(:,3,i))==0
            win_flag2 = 1;
        elseif max(B(:,4,i))==0
            win_flag2 = 1;
        elseif max(B(:,5,i))==0
            win_flag2 = 1;    
        end
    end
    
    %Check DIAGONALS and 4-CORNERS (both outer and inner)!
    if win_flag2 ==0
        
       %check backslash 
       m11 = B(1,1,i);
       m22 = B(2,2,i);
       m44 = B(4,4,i);
       m55 = B(5,5,i);
       t1 = [m11 m22 m44 m55];
       if max(t1) == 0
          win_flag2 = 1; 
       end
       
       %check forward slash
       m51 = B(5,1,i);
       m42 = B(4,2,i);
       m24 = B(2,4,i);
       m15 = B(1,5,i);
       t1 = [m15 m24 m42 m51];
       if max(t1) == 0
          win_flag2 = 1; 
       end
       
       %check OUTER 4-corners
       t1 = [m11 m15 m51 m55];
       if max(t1) == 0
           win_flag2 = 1;
       end
        
       %check INNER 4-corners
       t1 = [m22 m42 m24 m44];
       if max(t1) == 0
           win_flag2 = 1;
       end
       
    end
    
    if win_flag2 == 1
        num_winners = num_winners+1;
    end
    
    win_flag2 = 0;
    
end %ENDS LOOP OVER BOARDS
     
if num_winners > 0 
    win_flag = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Plays inner square bingo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B ct num_winners] = play_Inner_Square_Bingo(B,balls,ct)

win_flag = 0;   %initialize winning flag to end game
% ct: counter for number of balls drawn in regular bingo so far

while ( win_flag == 0 )

    ct = ct+1;        %updates counter when picking another ball
    ball = balls(ct); %draws a ball from already randomly drawn balls
    B = please_Dab_Bingo_Boards(B,ball); %dabs bingo board (sets elements to ZERO)
        
    %Checks for Regular Bingo
    [win_flag,num_winners] = please_Check_Inner_Square_For_Win(B); 
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Checks to see if any board has inner square completed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [win_flag,num_winners] = please_Check_Inner_Square_For_Win(B)

Nboards = length(B(1,1,:)); 
win_flag = 0;
num_winners = 0;
win_flag2 = 0;

for i=1:Nboards

    %check INNER SQUARE 
    m22 = B(2,2,i);
    m23 = B(2,3,i);
    m24 = B(2,4,i);
    m32 = B(3,2,i);
    m34 = B(3,4,i);
    m42 = B(4,2,i);
    m43 = B(4,3,i);
    m44 = B(4,4,i);
    t1 = [m22 m23 m24 m32 m34 m42 m43 m44];
    if max(t1) == 0
        win_flag2 = 1; 
    end
    
    if win_flag2 == 1
        num_winners = num_winners+1;
    end
    
    win_flag2 = 0;
    
end %ENDS LOOP OVER BOARDS

if num_winners > 0 
    win_flag = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Plays outer square bingo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B ct num_winners] = play_Outer_Square_Bingo(B,balls,ct)

first = 1;
win_flag = 0;   %initialize winning flag to end game
% ct: counter for number of balls drawn in regular bingo + inner square so far

while ( win_flag == 0 )
    
    if (first==1)
        
        %Checks for Outer Square Bingo BEFORE drawing another ball
        [win_flag,num_winners] = please_Check_Outer_Square_For_Win(B); 
        first = 0;
        
    else
        
        ct = ct+1;        %updates counter when picking another ball
        ball = balls(ct); %draws a ball from already randomly drawn balls
        B = please_Dab_Bingo_Boards(B,ball); %dabs bingo board (sets elements to ZERO)

        %Checks for Outer Square Bingo 
        [win_flag,num_winners] = please_Check_Outer_Square_For_Win(B); 
    end
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Checks if any boards have completed outer square bingo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [win_flag,num_winners] = please_Check_Outer_Square_For_Win(B)

Nboards = length(B(1,1,:)); 
win_flag = 0;
num_winners = 0;
win_flag2 = 0;

for i=1:Nboards

    %CHECK OUTER SQUARE 
    m11 = B(1,1,i);
    m12 = B(1,2,i);
    m13 = B(1,3,i);
    m14 = B(1,4,i);
    m15 = B(1,5,i);
    m51 = B(5,1,i);
    m52 = B(5,2,i);
    m53 = B(5,3,i);
    m54 = B(5,4,i);
    m55 = B(5,5,i);
    m21 = B(2,1,i);
    m31 = B(3,1,i);
    m41 = B(4,1,i);
    m25 = B(2,5,i);
    m35 = B(3,5,i);
    m45 = B(4,5,i);
    t1 = [m11 m12 m13 m14 m15 m51 m52 m53 m54 m55 m21 m31 m41 m25 m35 m45];
    if max(t1) == 0
        win_flag2 = 1; 
    end
    
    if win_flag2 == 1
        num_winners = num_winners+1;
    end
    
    win_flag2 = 0;
    
end %ENDS LOOP OVER BOARDS

if num_winners > 0 
    win_flag = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Plays cover-all bingo!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B ct num_winners] = play_Cover_All_Bingo(B,balls,ct)

first=1;
win_flag = 0;   %initialize winning flag to end game
% ct: counter for number of balls drawn in regular bingo + inner square so far

while ( win_flag == 0 )
    
    if (first == 1)
    
        %Checks for Cover All Bingo BEFORE drawing another ball
        [win_flag,num_winners] = please_Check_Cover_All_For_Win(B); 
        first = 0;    
        
    else
    
        ct = ct+1;        %updates counter when picking another ball
        ball = balls(ct); %draws a ball from already randomly drawn balls
        B = please_Dab_Bingo_Boards(B,ball); %dabs bingo board (sets elements to ZERO)
        
        %Checks for Cover All Bingo 
        [win_flag,num_winners] = please_Check_Cover_All_For_Win(B); 
   
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Checks to see if any board has completed cover all.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [win_flag,num_winners] = please_Check_Cover_All_For_Win(B)

Nboards = length(B(1,1,:)); 
win_flag = 0;
num_winners = 0;
win_flag2 = 0;

for i=1:Nboards

    %CHECK COVER ALL 
    if max(max(B(:,:,i))) == 0
        win_flag2 = 1; 
    end
    
    if win_flag2 == 1
        num_winners = num_winners+1;
    end
    
    win_flag2 = 0;
    
end %ENDS LOOP OVER BOARDS

if num_winners > 0 
    win_flag = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes averages for each set of monte carlo simulations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AvgMat = please_Compute_Averages(count,NgamesVec)

% count(j,i,k): 
% row,j:    game played for specific number of boards
% column,i: different number of board
% index, k: different number of simulations for convergence

AvgMat = zeros(length(NgamesVec),length(count(1,:,1)));

for k=1:length(NgamesVec)
   
    vec = sum(count(:,:,k)) / NgamesVec(k);
    AvgMat(k,:) = vec;
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Plots the number of games required to win all 4 bingo games
% (regular bingo, inenr square, outer square, cover all) in the most
% resolved Monte Carlo case (i.e., highest number of games played for each
% number of bingo cards)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Games_Required_To_Win(NgamesVec,BoardsVec,AvgMatReg,AvgMatIn,AvgMatOut,AvgMatCover)

%
% Plots all 4 different bingo games together 
%
figure(1);
sims = length(AvgMatCover(:,1));
plot(BoardsVec,AvgMatReg(sims,:),'r-'); hold on;
plot(BoardsVec,AvgMatIn(sims,:),'b-'); hold on;
plot(BoardsVec,AvgMatOut(sims,:),'k-'); hold on;
plot(BoardsVec,AvgMatCover(sims,:),'g-'); hold on;
plot(BoardsVec,AvgMatReg(sims,:),'r*'); hold on;
plot(BoardsVec,AvgMatIn(sims,:),'b*'); hold on;
plot(BoardsVec,AvgMatOut(sims,:),'k*'); hold on;
plot(BoardsVec,AvgMatCover(sims,:),'g*'); hold on;
strTitle = sprintf('Number of Draws Until You Win! (Averaging %d games per # of cards)', NgamesVec(end));
title(strTitle);
xlabel('Number of Bingo Cards');
ylabel('Average # of Draws to Win!');
legend('Regular Bingo','Inner Square','Outer Square','Cover All');

%
% Plots all 4 different bingo games separately
%
figure(2)
%
strLegend = sprintf('%d games averaged',NgamesVec(end));
%
subplot(2,2,1);
plot(BoardsVec,AvgMatReg(sims,:),'b-'); hold on;
plot(BoardsVec,AvgMatReg(sims,:),'r*'); hold on;
title('Regular Bingo');
xlabel('Number of Bingo Cards');
ylabel('Draws to Win');
legend(strLegend);
%
subplot(2,2,2);
plot(BoardsVec,AvgMatIn(sims,:),'b-'); hold on;
plot(BoardsVec,AvgMatIn(sims,:),'r*'); hold on;
title('Inner Square Bingo');
xlabel('Number of Bingo Cards');
ylabel('Draws to Win');
legend(strLegend);
%
subplot(2,2,3);
plot(BoardsVec,AvgMatOut(sims,:),'b-'); hold on;
plot(BoardsVec,AvgMatOut(sims,:),'r*'); hold on;
title('Outer Square Bingo');
xlabel('Number of Bingo Cards');
ylabel('Draws to Win');
legend(strLegend);
%
subplot(2,2,4);
plot(BoardsVec,AvgMatCover(sims,:),'b-'); hold on;
plot(BoardsVec,AvgMatCover(sims,:),'r*'); hold on;
title('Cover All Bingo');
xlabel('Number of Bingo Cards');
ylabel('Draws to Win');
legend(strLegend);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Plots the CONVERGENCE STUDY for how many games are required to
% get more accurate data for number of draws necessary to win at each type
% of bingo game.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Convergence_Study(BoardsVec,AvgMatReg,AvgMatIn,AvgMatOut,AvgMatCover)

%
% Plots CONVERGENCE STUDY for all 4 different bingo games separately 
%
figure(3);
subplot(2,2,1);
plot(BoardsVec,AvgMatReg(1,:),'m-'); hold on;
plot(BoardsVec,AvgMatReg(2,:),'b-'); hold on;
plot(BoardsVec,AvgMatReg(3,:),'k-'); hold on;
plot(BoardsVec,AvgMatReg(4,:),'g-'); hold on;
plot(BoardsVec,AvgMatReg(5,:),'r-','LineWidth',2); hold on;
plot(BoardsVec,AvgMatReg(1,:),'m-'); hold on;
plot(BoardsVec,AvgMatReg(2,:),'b-'); hold on;
plot(BoardsVec,AvgMatReg(3,:),'k-'); hold on;
plot(BoardsVec,AvgMatReg(4,:),'g-'); hold on;
plot(BoardsVec,AvgMatReg(1,:),'m*'); hold on;
plot(BoardsVec,AvgMatReg(2,:),'b*'); hold on;
plot(BoardsVec,AvgMatReg(3,:),'k*'); hold on;
plot(BoardsVec,AvgMatReg(4,:),'g*'); hold on;
plot(BoardsVec,AvgMatReg(5,:),'r*'); hold on;
title('Convergence Study: Reg. Bingo');
xlabel('Number of Bingo Cards');
ylabel('Average # of Draws to Win!');
legend('5','10','25','50','100');
%
subplot(2,2,2);
plot(BoardsVec,AvgMatIn(1,:),'m-'); hold on;
plot(BoardsVec,AvgMatIn(2,:),'b-'); hold on;
plot(BoardsVec,AvgMatIn(3,:),'k-'); hold on;
plot(BoardsVec,AvgMatIn(4,:),'g-'); hold on;
plot(BoardsVec,AvgMatIn(5,:),'r-','LineWidth',2); hold on;
plot(BoardsVec,AvgMatIn(1,:),'m-'); hold on;
plot(BoardsVec,AvgMatIn(2,:),'b-'); hold on;
plot(BoardsVec,AvgMatIn(3,:),'k-'); hold on;
plot(BoardsVec,AvgMatIn(4,:),'g-'); hold on;
plot(BoardsVec,AvgMatIn(1,:),'m*'); hold on;
plot(BoardsVec,AvgMatIn(2,:),'b*'); hold on;
plot(BoardsVec,AvgMatIn(3,:),'k*'); hold on;
plot(BoardsVec,AvgMatIn(4,:),'g*'); hold on;
plot(BoardsVec,AvgMatIn(5,:),'r*'); hold on;
title('Convergence Study: Inner Square');
xlabel('Number of Bingo Cards');
ylabel('Average # of Draws to Win!');
legend('5','10','25','50','100');
%
subplot(2,2,3);
plot(BoardsVec,AvgMatOut(1,:),'m-'); hold on;
plot(BoardsVec,AvgMatOut(2,:),'b-'); hold on;
plot(BoardsVec,AvgMatOut(3,:),'k-'); hold on;
plot(BoardsVec,AvgMatOut(4,:),'g-'); hold on;
plot(BoardsVec,AvgMatOut(5,:),'r-','LineWidth',2); hold on;
plot(BoardsVec,AvgMatOut(1,:),'m-'); hold on;
plot(BoardsVec,AvgMatOut(2,:),'b-'); hold on;
plot(BoardsVec,AvgMatOut(3,:),'k-'); hold on;
plot(BoardsVec,AvgMatOut(4,:),'g-'); hold on;
plot(BoardsVec,AvgMatOut(1,:),'m*'); hold on;
plot(BoardsVec,AvgMatOut(2,:),'b*'); hold on;
plot(BoardsVec,AvgMatOut(3,:),'k*'); hold on;
plot(BoardsVec,AvgMatOut(4,:),'g*'); hold on;
plot(BoardsVec,AvgMatOut(5,:),'r*'); hold on;
title('Convergence Study: Outer Square');
xlabel('Number of Bingo Cards');
ylabel('Average # of Draws to Win!');
legend('5','10','25','50','100');
%
subplot(2,2,4);
plot(BoardsVec,AvgMatCover(1,:),'m-'); hold on;
plot(BoardsVec,AvgMatCover(2,:),'b-'); hold on;
plot(BoardsVec,AvgMatCover(3,:),'k-'); hold on;
plot(BoardsVec,AvgMatCover(4,:),'g-'); hold on;
plot(BoardsVec,AvgMatCover(5,:),'r-','LineWidth',2); hold on;
plot(BoardsVec,AvgMatCover(1,:),'m-'); hold on;
plot(BoardsVec,AvgMatCover(2,:),'b-'); hold on;
plot(BoardsVec,AvgMatCover(3,:),'k-'); hold on;
plot(BoardsVec,AvgMatCover(4,:),'g-'); hold on;
plot(BoardsVec,AvgMatCover(1,:),'m*'); hold on;
plot(BoardsVec,AvgMatCover(2,:),'b*'); hold on;
plot(BoardsVec,AvgMatCover(3,:),'k*'); hold on;
plot(BoardsVec,AvgMatCover(4,:),'g*'); hold on;
plot(BoardsVec,AvgMatCover(5,:),'r*'); hold on;
title('Convergence Study: Cover All');
xlabel('Number of Bingo Cards');
ylabel('Average # of Draws to Win!');
legend('5','10','25','50','100');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints simulation information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_info()


fprintf('\n\n**********************************************************************************\n');
fprintf('\nAuthor: Nicholas A. Battista \n');
fprintf('University: UNC-CH \n');
fprintf('First Created: January 19, 2015 \n');
fprintf('Last Revision: May 5, 2015 \n\n');
fprintf('**********************************************************************************\n\n');

fprintf('This code plays the game of Bingo in a Monte Carlo setting, averaging\n');
fprintf('   N simulated games of Bingo for a specified number of bingo boards \n\n');

fprintf('It tries to answer the following questions: \n');
fprintf('    1. What is the expected number of draws to win regular bingo\n');
fprintf('       if there are N bingo boards in play? \n');
fprintf('    2. What is the convergence rate of expected number of draws \n');
fprintf('       vs. # of bingo boards for a certain number of simulated bingo games? \n');
fprintf('    3. What about inner square, outer square, and cover all? \n\n');

fprintf('Initially it will play {5,10,25,50,100} games of bingo with each number\n');
fprintf('    of bingo boards. It varies the number of boards between {1,5,10,15,...,75}.\n\n');

fprintf('NOTE: Monte Carlo simulation will take roughly 5 minutes. Please sit back and enjoy!\n'); 

fprintf('\n**********************************************************************************\n');
