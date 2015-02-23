clear all

nbrOfAgents = 1500;
nbrOfStrategies = 2; %Number of Strategies
slutM=4; %Largest memory
nbrOfTimeSteps = 12000; %Time steps of one game
nbrOfRuns = 1; %Number of runs with each setup
epsilon=0.1; %Tradingfee
lambda=1;%A strange parameter, something with liquidity
partProducer=0.667; %Ratio of producer/nbrOfAgents

Prod=floor(partProducer*nbrOfAgents); %Number of Producers


%Pre-allocating some measured quantities
sigma=ones(nbrOfRuns,slutM);
price=ones(nbrOfTimeSteps+1,1)*1000;
priceIncr=zeros(nbrOfTimeSteps,1);
stratCompare=zeros(nbrOfTimeSteps,2);

for memory = slutM:slutM
            P=2^memory;
    % Run simulation
    %according to "Black-Scholes to MG" GCMG depends on an alpha_c for both 
    %Producers and speculators.
    
    for iRun = nbrOfRuns:nbrOfRuns
        
        agents = CanonicalAgents(nbrOfAgents,memory,nbrOfStrategies); %Creates a number
        %of agents which consists of 3 indices with (whichAgents,Matrix with strategy outcomes (2^m,s))
        
        points = (epsilon+1)*ones(nbrOfAgents, nbrOfStrategies); %A separate matrix which
        %corresponds to the agents' matrix
        
        storeMinorityGroup = zeros(nbrOfTimeSteps,1);
        storeAttendance = zeros(nbrOfTimeSteps,1); %Resetting
        
        for timeStep = 1:nbrOfTimeSteps
            
            actions = zeros(nbrOfAgents,1);%Resetting
            
            histIndex = randi(2^memory,1,1);%Rand one of the 2^memory stra
            %tegies
            
            
            % Which will go
            %Selects the choice perscribed by the at this point most
            %succesful strategy. If there is more than one at the top
            %a coin toss decides the outcome.
            %With a possibility to abstained if the threshold of a
            %stratey is not met.
            
            %This could be modified with one loop and different epsilon
            %instead.
            for i=1:Prod
                
                actions(i)=agents(i,histIndex,1);
            end
            %Perhaps no coin toss if two strategies are equally good
            %instead always choose the same one.
            %Would make analysis easier and irrational agents (such)
            %as ours are likely to act on a hunch.
            for i=Prod+1:nbrOfAgents
                if max(points(i,:))>0
                    
                    [Value,stratIndex] = max(points(i,:),[],2);
                    maxIndex = find(points(i,:)==Value);
                    
                    if(length(maxIndex)>1)
                        stratIndex = maxIndex(floor(rand*length(maxIndex))+1);
                    end
                    
                    actions(i)=agents(i,histIndex,stratIndex);
                    
                else
                    actions(i)=0;
                end
                
            end
            
            
            minority=-sign(sum(actions));
            
            %           %Prisättning med log
            %             price(timeStep+1)=price(timeStep)*exp(sum(actions)/lambda);
            %             priceIncr(timeStep)=exp(sum(actions)/lambda);
            
            
            %Difference of sum actions determines the price increase.
            price(timeStep+1)=price(timeStep)*exp(sum(actions)/sum(abs(actions)));
            %priceIncr(timeStep)=sum(actions);
            
            
            for i=1:nbrOfAgents %Update score
                
                for j=1:nbrOfStrategies %for each strategy
                    points(i,j) = points(i,j) - agents(i,histIndex,j)*sum(actions)/P-epsilon/P;
                    %points(i,j) = points(i,j) + agents(i,histIndex,j)*minority;
                    
                end
            end
            stratCompare(timeStep,:)=points(nbrOfAgents,:);            
            
            storeAttendance(timeStep) = sum(abs(actions)); %Number of agents
            %that doesn't abstain.
            
                     
        end %timestep
        %Store the measure of standard deviation.
        sigma(iRun,memory)=std(storeMinorityGroup);
        
        
    end % End of Runs
        
end     %End memory

%% Plottar lite goe grejer

figure(1)
clf
subplot(2,1,1)
title('Participants')
hold on
%plot(storeAttendance(1:end))
%axis([0 nbrOfTimeSteps 0 nbrOfAgents+5])
activeSpec=storeAttendance-floor(partProducer*nbrOfAgents);

[val ind]=min(activeSpec);
width=1000;
start=ind-width/2;
stopp=ind+width/2;
%plot(activeSpec(start:stopp))
plot(activeSpec)
subplot(2,1,2)
title('Price (log)')
% plot(price(start:stopp))
plot(price)

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
 

plot(logReturn2)

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
x=1:maxLag;

coeff= coeffvalues(fit((1:maxLag)',corr,'Power1'))
theorCorr=coeff(1)*x.^coeff(2);

figure(3)
clf
plot(theorCorr,'--')
hold on
plot(corr,'o')
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
plot(logReturn2(4000:end))

subplot(2,1,2)
hold on
title('Fördelning av prisökning')
hist(logReturn2(4000:end),20)
