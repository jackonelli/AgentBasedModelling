%% Main-fil för simulering
cd /Users/jakoblindqvist/GitHub/AgentBasedModeling

%%

%************
% PARAMETERS
%************
nbrOfAgents = 300;
nbrOfStrategies = 2;
nbrOfTimeSteps = 20000;
nbrOfRuns = 1;
epsilon=0.01;
memory=5;

a=1;

partProducer=0.05;    %Ratio of producer/nbrOfAgents
%Alpha_c


%****************
% RUN SIMULATION
%****************
for memory=4:4
    P=2^memory;
    %Secondary parameters
    Prod=floor(partProducer*nbrOfAgents); %Number of Producers
    
    %*********************************
    %PRE-ALLOCATING STORED QUANTITIES
    %*********************************
    
    activeSpec = zeros(nbrOfTimeSteps,nbrOfRuns);
    actions = zeros(nbrOfAgents-Prod,1);
    A = zeros(nbrOfTimeSteps,nbrOfRuns);
    bestStrat=zeros(nbrOfAgents-Prod,nbrOfTimeSteps);
    clear agents
    
    AProd=[-2 -2 -2 -2 -1 -1 0 0 0 0 1 1 2 2 2 4];
    %Deterministic Producers for m=4
    
    %AProd=[(-P/2:1) (1:P/2)];
    %AProd(1)=max(AProd)+2;
    
    for iRun=1:nbrOfRuns
        
        %AProd=sum(ones(Prod,P) - floor(rand(Prod,P)*2)*2);
        %store_Prod(iRun,:)=AProd;
        %Random contribution to the A(t) variable by the producers Correct stat.
        %distr. but not code efficient
        storeHist=zeros(nbrOfTimeSteps,1);
        
        %Resetting dynamic values
        agents = CanonicalAgents(nbrOfAgents-Prod,P,nbrOfStrategies);
        %storeAgree(iRun,1)=sum(sum((squeeze(agents(:,5,:)*AProd(5)>0)),2)==2)
        points = ones(nbrOfAgents-Prod,1, nbrOfStrategies);
        %A separate matrix which corresponds to the agents' matrix
        
        for timeStep = 1:nbrOfTimeSteps
            %Rand one of the P=2^memory histories
            %Controlling when the extreme history comes (at index P) in
            %order to control bubbles forming
%             if (timeStep>8000 && max(bestStrat(:,timeStep-1))>40);
%                 histIndex = randi(P-1,1,1);
%             else
                histIndex = randi(P,1,1);
%             end
%             
            if(histIndex==P);
                storeHist(timeStep)=1;
            end
%             
            %**********************
            % WHICH AGENTS WILL GO
            %**********************
            
            %Selects the choice perscribed by the at this point most
            %succesful strategy. If there is more than one at the top
            %a coin toss decides the outcome.
            %With a possibility to abstain if the threshold of a
            %strategy is not met.
            
            [Value,stratIndex] = max(points,[],3);
            %Pick most successfull strategy and it's corresp. point.
            tempStrat=agents(:,histIndex,:);
            %Intermediate step because Matlab doeesn't do 3D-matrices well
            
            % actions=tempStrat(stratIndex).*Value>0;%
            %Jakobs finlösning: fungerar ej
            
            %             ind=find(points(:,1)==points(:,2)); %Samuel hävdar att find inte behövs
            %             temp2=tempStrat(:,1).*(stratIndex==1)+tempStrat(:,2).*(stratIndex==2);
            %             temp2(ind)=ones(length(ind),1)-floor(rand(length(ind),1)*2)*2;
            
            ind=(points(:,1,1)==points(:,1,2));
            temp2=tempStrat(:,1,1).*(stratIndex==1)+tempStrat(:,1,2).*(stratIndex==2);
            temp2(ind)=ones(sum(ind),1)-floor(rand(sum(ind),1)*2)*2;
            
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
            
            
            activeSpec(timeStep,iRun)=sum(abs(actions));
            A(timeStep,iRun)=(sum(actions)+AProd(histIndex));%/(activeSpec(timeStep,iRun)+Prod);
            %**************
            % UPDATE SCORE
            %**************
            
            points = points(:,1,:) - tempStrat*((sum(actions)+AProd(histIndex)))-epsilon;
            
            bestStrat(:,timeStep)=Value;
            
        end % end timestep
        
    end
    
end

%%
figure(1)
clf
plot(activeSpec,'-')
hold on
plot(storeHist*40,'k.')
%
figure(2)
clf
for i=1:nbrOfAgents-Prod
    plot(bestStrat(i,:))
    hold on
end
plot(linspace(0,nbrOfTimeSteps),-epsilon*linspace(0,nbrOfTimeSteps),'r--','LineWidth',1)
plot(storeHist*80,'.k')
set(gca,'FontSize',18)
xlabel('Tidssteg')
ylabel('Poäng')

%%
epsPhase=7000;
bub=find(activeSpec(epsPhase:end)>40)+epsPhase;
B1=diff(bub)>40;
B1(1)=1;
figure(3)
clf
plot(bub,'r+')
hold on
plot(bub(B1),'b*')

%% *******************
%PRICE AND LOGRETURNS
%  *******************
epsPhase=10000;
stopp=nbrOfTimeSteps;

price=1000*ones(nbrOfTimeSteps+1,1);
lambda=4E3;
for timeStep=1:nbrOfTimeSteps
    price(timeStep+1)=price(timeStep)*exp(A(timeStep)/lambda);
end

%Calculate log Returns
dt=1;
%logReturn=log((price(dt+1:end)-price(1:end-dt))./price(1:end-dt));
%New logreturn that doesn't work

%Return=log(price(dt+1:end)./price(1:end-dt));
%NBins=100;
Return=A(epsPhase:stopp,1)/lambda;

distrFit(Return,10000,'gauss1');


figure(2)
clf
plot(smooth(price(epsPhase:end),100))


%% ****************
%  Autocorrelation
%  ****************
clear corr
maxLag=1000;
epsPhase=10000;
corr=autocorr(abs(A(epsPhase:end)),maxLag);
corr=corr(2:end);
x=(1:maxLag)';


[powerCorr,gof]=fit((1:maxLag)',corr,'Power1');
coeff=coeffvalues(powerCorr);
%powerCoeff= coeffvalues(fit((1:maxLag)',corr,'Power1'));
%expCoeff = coeffvalues(fit((1:maxLag)',corr,'Exp1'));


%powerCorr=powerCoeff(1)*x.^powerCoeff(2);
%expCorr=expCoeff(1)*exp(expCoeff(2)*x);

figure(3)
clf
plot(powerCorr,x,corr)
set(gca,'FontSize',16)
title('Autocorrelated logReturn')
legend(strcat('b=-',num2str(coeff(2))))
%legend('Autocorr',strcat('ax^{b}, b=',num2str(powerCoeff(2))),strcat('ae^{bx}, b=',num2str(expCoeff(2))))
%% Autocorrelation for different alpha



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


%% *********
%  PRICE GIF
%  *********
dt=10;
I=nbrOfTimeSteps/dt;
tStart=12000;

clf
for i=1:I
    figure(1)
    r=tStart+(i*dt:(i+1)*dt);
    plot(r,price(r),'b-');
    xlim([0 nbrOfTimeSteps-tStart])
    ylim([500 1300])
    hold on
    pause(0.01)
end

