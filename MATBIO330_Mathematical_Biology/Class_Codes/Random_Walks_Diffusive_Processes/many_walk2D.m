function many_walk2D(N, M);

%Simulates many random walkers that have an equally likely chance of moving
%left, right, up or down.
%N is the number of time steps
%M is the number of walkers
clf;

x=zeros(1,N); %initialize the vector to hold the x-positions
y=zeros(1,N); %initialize the vector to hold the y-positions
x(1,1)=0;   %set initial x-position of walker
y(1,1)=0;   %set initial y-position of walker

%set parameters here
stepsize = 1.0;   %stepsize of walker
bins= N/80;   %size of bins
binmax = N/2;   %maximum bin
binnum = round(binmax/bins);    %number of bins

rad = zeros(1,M);  %final distance of walker
nn = 0:bins:binmax;  %middle value gives bin size for histogram
size1=size(nn);
clf;

for j=1:M, %loop through all walkers
    for i=2:N, %loop through all of the time steps
        r = rand; %pick a random number to determine up/down or right/left
        if r>0.5, %go right/left
            r1 = rand; %pick a random number to determine right or left
            if r1>0.5, %go right
                y(1,i)=y(1,i-1);
                x(1,i)=x(1,i-1)+stepsize;    
            else %go left
                y(1,i)=y(1,i-1);
                x(1,i)=x(1,i-1)-stepsize;
            end
        else %go up/down
            r2=rand; %pick a random number to go up or down
            if r2>0.5, %go up
                x(1,i)=x(1,i-1);
                y(1,i)=y(1,i-1)+stepsize;    
            else %go down
                x(1,i)=x(1,i-1);
                y(1,i)=y(1,i-1)-stepsize;
            end
        end
    end

    rad(1,j)=sqrt(x(1,N)^2+y(1,N)^2); %This holds the final distance from the starting point of each walker
    
    %use these to set the axes
    xmax = max(abs(x))+1;
    ymax = max(abs(y))+1;

    figure(1)
    hold on
    plot(x,y,'r');    % plot the path
    axis([-N N -N N]);
    axis square
    title(['Two-dimensional ' int2str(N) '-step random walk']);
    zoom on
    
end

%calculate theoretical distribution
x1=zeros(1,100);
pr=zeros(1,100);
dx=nn(size1(2))/100;

for k=1:100,
    x1(1,k)=(k-1)*dx;
    pr(1,k)= (2*x1(1,k)/((stepsize^2)*N))*exp(-x1(1,k)*x1(1,k)/((stepsize^2)*N));
end

pr=pr.*(bins*M);   %rescale theoretical distribution

figure(2)
%plot histogram and theoretical distribution
hist(rad,nn)
hold on
plot(x1,pr, 'r')
xlabel('Distance')
ylabel('Number of walkers')