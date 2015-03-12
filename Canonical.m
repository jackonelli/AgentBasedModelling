%% Main-fil för simulering
% /Users/jakoblindqvist/GitHub/AgentBasedModeling

%%
clear all

%************
% PARAMETERS
%************
nbrOfAgents = 501;
nbrOfStrategies = 2;
nbrOfTimeSteps = 20000;
nbrOfRuns = 1;
epsilon=0.1;          %Tradingfee
partProducer=0.2;    %Ratio of producer/nbrOfAgents
memory=4;

%Secondary parameters
Prod=floor(partProducer*nbrOfAgents); %Number of Producers
AProd=sum(ones(Prod,2^memory) - floor(rand(Prod,2^memory)*2)*2);
%Random contribution to the A(t) variable by the producers Correct stat. 
%distr. but not code efficient

P=2^memory;

%*********************************
%PRE-ALLOCATING STORED QUANTITIES
%*********************************

price=ones(nbrOfTimeSteps+1,1)*1000;
priceIncr=zeros(nbrOfTimeSteps,1);
stratCompare=zeros(nbrOfTimeSteps,2);
activeSpec = zeros(nbrOfTimeSteps,nbrOfRuns);
A = zeros(nbrOfTimeSteps,nbrOfRuns);
actions = zeros(nbrOfAgents-Prod,1);
bestStrat=zeros(nbrOfAgents-Prod,nbrOfTimeSteps);

%****************
% RUN SIMULATION
%****************

for iRun=1:nbrOfRuns
    
    %Resetting dynamic values
    agents = CanonicalAgents(nbrOfAgents-Prod,memory,nbrOfStrategies);
    
    points = (epsilon)*ones(nbrOfAgents-Prod, nbrOfStrategies);
    %A separate matrix which corresponds to the agents' matrix
    
    history=zeros(1,nbrOfTimeSteps+memory-1);
    history(1:memory)=randi(2,1,memory)-1;
    %Initiate history vector and create a random first history
    
    for timeStep = 1:nbrOfTimeSteps
        
        histIndex=bi2de(history(1,timeStep:timeStep+memory-1))+1;
        
        
        %Rand one of the 2^memory histories
        
        %**********************
        % WHICH AGENTS WILL GO
        %**********************
        
        %Selects the choice perscribed by the at this point most
        %succesful strategy. If there is more than one at the top
        %a coin toss decides the outcome.
        %With a possibility to abstain if the threshold of a
        %strategy is not met.
        
        [Value,stratIndex] = max(points,[],2);
        % Pick most successfull strategy and its corresp. point.
        tempStrat(:,:)=agents(:,histIndex,:);

        %Intermediate step because Matlab doeesnt do 3D-matrices well
        
        % actions=tempStrat(stratIndex).*Value>0;%
        %Jakobs finlösning: fungerar ej
        
        
        
        ind=find(points(:,1)==points(:,2));
        temp2=tempStrat(:,1).*(stratIndex==1)+tempStrat(:,2).*(stratIndex==2);
        temp2(ind)=ones(length(ind),1)-floor(rand(length(ind),1)*2)*2;
        %Samuels fullösning: FUNGERAR.
        
        actions=temp2.*(Value>0);
        %Stores the actual choice of all speculators.
        
        %(This could be modified with one loop and different epsilon
        %instead.
        %Perhaps no coin toss if two strategies are equally good
        %instead always choose the same one.
        %Would make analysis easier and irrational agents (such)
        %as ours are likely to act on a hunch.)
        
        %**************
        %STORE RESULTS
        %**************
        
        A(timeStep,iRun)=sum(actions)+AProd(histIndex);
        activeSpec(timeStep,iRun)=sum(abs(actions));
        
        %**************
        % UPDATE SCORE
        %**************
        
        points = points - tempStrat*A(timeStep)/P-epsilon/P;
        bestStrat(:,timeStep)=Value;
        
        if A(timeStep,iRun)==0
        history(timeStep)=randi(2)-1;
        else
        history(timeStep)=1==-sign(A(timeStep,iRun));
        end
        
    end % end timestep
end

%% **************
%  BEST STRATEGY
%  **************

figure(3)
clf
for a=1:(nbrOfAgents-Prod)
    hold on
    plot(bestStrat(a,:))
end

% *******************
%PRICE AND LOGRETURNS;
%  *******************

%Calculate price
price=1000*ones(nbrOfTimeSteps+1,1);
lambda=1E5;
for timeStep=1:nbrOfTimeSteps
    price(timeStep+1)=price(timeStep)*exp(A(timeStep)/lambda);
end

%Calculate log Returns
dt=1;
logReturn=log(price(dt+1:end)./price(1:end-dt));

%plot
% figure(1)
% clf
% hold on
% title('Price')
% xlabel('Time')
% plot(price)

figure(2)
subplot(2,1,1)
plot(activeSpec(1:end))
title('')
%hist(logReturn(1:end),500)

subplot(2,1,2)
plot(logReturn(1:end))
title('log Returns')
%% ****************
%  Autocorrelation
%  ****************
clear corr
maxLag=300;
corr=autocorr(abs(logReturn(2000:end)),maxLag);
corr=corr(2:end);

coeff= coeffvalues(fit((1:maxLag)',corr,'Power1'))
x=linspace(1,maxLag);
theorCorr=coeff(1)*x.^coeff(2);

figure(3)
clf
hold on
plot(theorCorr,'--')
plot(corr,'o')
title('Autocorrelated logReturn')
% Corr is the vector of auto correlations.
%loglog(corr,'o')



