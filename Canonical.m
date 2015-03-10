%% Main-fil för simulering
% /Users/jakoblindqvist/GitHub/AgentBasedModeling

%%
tic
nbrOfAgents = 501;
nbrOfStrategies = 2;    %Number of Strategies
%Largest memory
nbrOfTimeSteps = 10000; %Time steps of one game
nbrOfRuns = 4;          %Number of runs with each setup
epsilon=0.1;          %Tradingfee
lambda=1;               %A strange parameter, something with liquidity
partProducer=0.5;     %Ratio of producer/nbrOfAgents

Prod=floor(partProducer*nbrOfAgents); %Number of Producers
memory=4;

%Pre-allocating some measured quantities
price=ones(nbrOfTimeSteps+1,1)*1000;
priceIncr=zeros(nbrOfTimeSteps,1);
stratCompare=zeros(nbrOfTimeSteps,2);



% Run simulation
%according to "Black-Scholes to MG" GCMG depends on an alpha_c for both
%Producers and speculators.
P=2^memory;
AProd=sum(ones(Prod,2^memory) - floor(rand(Prod,2^memory)*2)*2);
%Correct stat. distr. but not code efficient
activeSpec = zeros(nbrOfTimeSteps,nbrOfRuns); %Resetting

A = zeros(nbrOfTimeSteps,1);
actions = zeros(nbrOfAgents-Prod,1);%Resetting

for iRun=1:nbrOfRuns
    
    agents = CanonicalAgents(nbrOfAgents-Prod,memory,nbrOfStrategies); %Creates a number
    %of agents which consists of 3 indices with (whichAgents,Matrix with strategy outcomes (2^m,s))
    points = (epsilon)*ones(nbrOfAgents-Prod, nbrOfStrategies); %A separate matrix which
    %corresponds to the agents' matrix
    clear temp1; clear tempStrat; clear temp;
    %rr = randi(2,nbrOfAgents-Prod,nbrOfTimeSteps);
    for timeStep = 1:nbrOfTimeSteps
        
        
        
        histIndex = randi(2^memory,1,1);%Rand one of the 2^memory histories
        
        %**********************
        % WHICH AGENTS WILL GO
        %**********************
        
        %Selects the choice perscribed by the at this point most
        %succesful strategy. If there is more than one at the top
        %a coin toss decides the outcome.
        %With a possibility to abstained if the threshold of a
        %strategy is not met.
        
        
        [Value,stratIndex] = max(points,[],2); %Pick most success-
        %                 %full strategy and it's corresp. point.
        tempStrat(:,:)=agents(:,histIndex,:); %Intermediate step
        %                 %because Matlab doeesn't do 3D-matrices well
        
        
        
        % actions=tempStrat(stratIndex).*Value>0;%
        %Jakobs finlösning: fungerar ej
        
        
        
        ind=find(points(:,1)==points(:,2));
        temp2=tempStrat(:,1).*(stratIndex==1)+tempStrat(:,2).*(stratIndex==2);
        temp2(ind)=ones(length(ind),1)-floor(rand(length(ind),1)*2)*2;
        %notEqual = points(:,1)~=points(:,2);
        %r = randi(2,nbrOfAgents-Prod,1);
        %r= rr(:,timeStep);
        %temp2 = (tempStrat(:,1).*(stratIndex==1) + tempStrat(:,2).*(stratIndex==2)).*notEqual +...
           % (ones(size(points,1),1)-notEqual).*(tempStrat(:,1).*(r==1) + tempStrat(:,2).*(r==2));
        
        actions=temp2.*(Value>0);
        %Samuels fullösning: FUNGERAR.
        
        %Stores the actual
        %choice of all speculators.
        
        
        %(This could be modified with one loop and different epsilon
        %instead.
        
        %Perhaps no coin toss if two strategies are equally good
        %instead always choose the same one.
        %Would make analysis easier and irrational agents (such)
        %as ours are likely to act on a hunch.)
        
        %**************
        %STORE RESULTS
        %**************
        
        A(timeStep)=sum(actions)+AProd(histIndex);
        activeSpec(timeStep,iRun)=sum(abs(actions));
        %*************
        % UPDATE SCORE
        %**************
        
        points = points - tempStrat*A(timeStep)/P-epsilon/P;
        %points(i,j) = points(i,j) + agents(i,histIndex,j)*minority;
        
        
        
    end % end timestep
    
end %end nbrOfRuns
toc
%% Tidtagning
tidsst=[100 500 1000 5000 10000 20000 40000];
tid=[0.043 0.14 0.26 1.15 2.3 5.57 9.11];
figure(1)
clf
plot(tidsst,tid)

%%
figure(1)
clf
for i=1:4
    subplot(2,2,i)
    plot(activeSpec(:,i))
end

%% Plottar lite goe grejer

figure(1)
clf
subplot(2,1,1)
title('Participants')
hold on
%plot(storeAttendance(1:end))
%axis([0 nbrOfTimeSteps 0 nbrOfAgents+5])
activeSpec=storeAttendance-floor(partProducer*nbrOfAgents);

[val,ind]=min(activeSpec);
width=1000;
start=ind-width/2;
stopp=ind+width/2;
%plot(activeSpec(start:stopp))

plot(activeSpec)

subplot(2,1,2)
hold on
title('Price')
% plot(price(start:stopp))
plot(price)

%% Price

price(timeStep+1)=price(timeStep)*exp(sum(actions)/sum(abs(actions)));

%% Logreturn
dt=1;
clear logReturn
clear logReturn2
figure(2)
clf


for i=1+dt:nbrOfTimeSteps
    
    logReturn(i)=log(price(i)/price(i-dt));
    
end
logReturn2=logReturn(1,(dt+1):end);


plot(logReturn2(2000:end))

%
% for i=1+dt:nbrOfTimeSteps
% avkast(i)=(price(i)-price(i-dt))/(price(i));
% end
% avkast2=avkast(1,2:end)';
%
%% Autocorrelation
clear corr
maxLag=300;
corr=autocorr(abs(logReturn2(2000:end)),maxLag);
corr=corr(2:end)';


coeff= coeffvalues(fit((1:maxLag)',corr,'Power1'))
theorCorr=coeff(1)*x.^coeff(2);

figure(3)
clf
hold on
plot(theorCorr,'--')
plot(corr,'o')
title('Autocorrelated logReturn')
% Corr is the vector of auto correlations.
%loglog(corr,'o')

%%
figure(4)
clf
k=1.1;
subplot(2,1,1)
title('logReturn')
%axis([0,length(logReturn2),-k*max(logReturn2),k*max(logReturn2)])
hold on
plot(logReturn2(2:end))

subplot(2,1,2)
hold on
title('Fördelning av prisökning')
hist(logReturn2(2:end),100)



