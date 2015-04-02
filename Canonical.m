%% Main-fil för simulering
% /Users/jakoblindqvist/GitHub/AgentBasedModeling

%%


%************
% PARAMETERS
%************
nbrOfAgents = 500;
nbrOfStrategies = 2;
nbrOfTimeSteps = 100000;
nbrOfRuns = 1;
epsilon=0.01;          %Tradingfee
% partProducer=0.4;    %Ratio of producer/nbrOfAgents
memory=4;

partProducer=0.3;


%Alpha_c
alpa=0;
sigmaScaled=0;
phase=5000;


%*********************************
%PRE-ALLOCATING STORED QUANTITIES
%*********************************

%price=ones(nbrOfTimeSteps+1,1)*1000;
%priceIncr=zeros(nbrOfTimeSteps,1);


%****************
% RUN SIMULATION
%****************
for memory=4:4
    activeSpec = zeros(nbrOfTimeSteps,nbrOfRuns);
    actions = zeros(nbrOfAgents-Prod,1);
    %bestStrat=zeros(nbrOfAgents-Prod,nbrOfTimeSteps);
    A = zeros(nbrOfTimeSteps,nbrOfRuns);
    clear agents
    
    
    P=2^memory;
    %Secondary parameters
    Prod=floor(partProducer*nbrOfAgents); %Number of Producers
    
    %AProd=sum(ones(Prod,P) - floor(rand(Prod,P)*2)*2);
    %Random contribution to the A(t) variable by the producers Correct stat.
    %distr. but not code efficient
    
    AProd=[-7 -4 -4 -3 -2 -1 0 0 0 0 1 2 3 4 4 7];
    %Deterministic Producers for m=4
    
    for iRun=1:nbrOfRuns
        
        %Resetting dynamic values
        agents = CanonicalAgents(nbrOfAgents-Prod,memory,nbrOfStrategies);
        %storeAgree(iRun,1)=sum(sum((squeeze(agents(:,5,:)*AProd(5)>0)),2)==2)
        points = zeros(nbrOfAgents-Prod, nbrOfStrategies);
        %A separate matrix which corresponds to the agents' matrix
        
        for timeStep = 1:nbrOfTimeSteps
            
            histIndex = randi(P,1,1);
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
            %Pick most successfull strategy and it's corresp. point.
            tempStrat=squeeze(agents(:,histIndex,:));
            %Intermediate step because Matlab doeesn't do 3D-matrices well
            
            % actions=tempStrat(stratIndex).*Value>0;%
            %Jakobs finlösning: fungerar ej
            
            ind=find(points(:,1)==points(:,2));
            temp2=tempStrat(:,1).*(stratIndex==1)+tempStrat(:,2).*(stratIndex==2);
            temp2(ind)=ones(length(ind),1)-floor(rand(length(ind),1)*2)*2;
            %Samuels fullösning: FUNGERAR.
            
            actions=temp2.*(Value>0);
            %Stores the actual choice of all speculators.ind
            
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
            
            points = points - tempStrat*A(timeStep,iRun)-epsilon;
            %bestStrat(:,timeStep)=Value;
            
        end % end timestep
        
    end
    
    
    
end

%% *******************
%PRICE AND LOGRETURNS
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

start=1;
stopp=nbrOfTimeSteps;

figure(1)
clf
hold on
title('LogReturns')
xlabel('Time')
plot(A(start:stopp))

figure(2)
clf
subplot(2,1,1)
plot(activeSpec(start:stopp,1))
title('Active Speculators')
%hist(logReturn(1:end),500)

subplot(2,1,2)
plot(price(start:stopp))
title('Price')
%%
figure(3)
clf
hist(price(1E4,end),10000)
%% ***************************
%ALPHA_C scaled with N_spec.
%***************************

%     volA(memory)=mean(std(abs(A(phase:end))))^2;


%************************************
% ALPHA_C WITH MEAN AFTER PHASE SHIFT
%************************************
%This is probably shit because there is conflict between stat. values
%and deterministic ones.
%     AStab=A(phase:end,:);
%     specStab=mean(activeSpec(phase:end,:));
%     alpa(memory)=P./mean(specStab)
%     sigma(memory)=mean(std(AStab).^2./specStab);
%



%*******************************
% ALPHA_C WITH FLOATING MEAN/STD
%*******************************
    timeInt=600;
    A=A(6001:end);
    nAct=zeros(length(A)/timeInt,1);
    stdA=zeros(size(nAct));

    for i=0:length(A)/timeInt-1

        ind1=i*timeInt+1;
        ind2=(i+1)*timeInt;

        nAct(i+1)=mean(activeSpec(ind1:ind2))+Prod;
        stdA(i+1)=std(abs(A(ind1:ind2)));

        %Perhaps avoid for-loop here with a cunning mean function

    end

    alpa=[alpa; 2^memory./nAct/nbrOfAgents^(-0.65)];
    sigmaScaled=[sigmaScaled; stdA.^2./nAct];

% Alpha_c, Weighted with activeSpec.

figure(1)
hold on
loglog(alpa,sigmaScaled,'.')
%loglog(alpa,volScaled,'.')
xlabel('\alpha=2^m/n_a^s','FontSize',18)
ylabel('\sigma^2/n_a^s','FontSize',18)

%% Colour plot with all best strategies
figure(1)
clf
for i=1:nbrOfAgents-Prod
    plot(bestStrat(i,:))
    hold on
end

%%
phase=2;
figure(2)
clf
plot(activeSpec(1:end))
minim=min(activeSpec(phase:end))
maxim=max(activeSpec(phase:end))
%medel=mean(activeSpec(phase+5000:phase+7000))
%P/medel


%% *******************
%PRICE AND LOGRETURNS
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

start=1;
stopp=30000;

figure(1)
clf
hold on
title('LogReturns')
xlabel('Time')
plot(A(start:stopp))

figure(2)
clf
subplot(2,1,1)
plot(activeSpec(start:stopp,1))
title('Active Speculators')
%hist(logReturn(1:end),500)

subplot(2,1,2)
plot(price(start:stopp))
title('Price')
%% ****************
%  Autocorrelation
%  ****************
clear corr
maxLag=300;
corr=autocorr(abs(A(start:end)),maxLag);
corr=corr(2:end);

coeff= coeffvalues(fit((1:maxLag)',corr,'Power1'))
x=linspace(1,maxLag);
theorCorr=coeff(1)*x.^coeff(2);

figure(3)
clf
hold on
plot(x,theorCorr,'--')
plot(corr,'o')
title('Autocorrelated logReturn')
% Corr is the vector of auto correlations.
%loglog(corr,'o')



