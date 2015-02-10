clear all

nbrOfAgents = 101;
s = 2; %Number of Strategies
slutM=2; %Largest memory
nbrOfTimeSteps = 1500; %Time steps of one game
nbrOfRuns = 1; %Number of runs with each setup
epsilon=0.2; %Tradingfee
lambda=1.5;%A strange parameter, something with liquidity

%Pre-allocating of some measured quantities
sigma=ones(nbrOfRuns,slutM);
price=ones(nbrOfTimeSteps+1,1)*1000;
priceIncr=zeros(nbrOfTimeSteps,1);

for memory = slutM:slutM
    
    % Run simulation
    
    for iRun = nbrOfRuns:nbrOfRuns
        
        agents = CanonicalAgents(nbrOfAgents,memory,s); %Creates a number
        %of agents which consists of 3 indices with (whichAgents,Matrix with strategy outcomes (2^m,s))
        
        points = (epsilon+1)*ones(nbrOfAgents, s); %A separate matrix which
        %corresponds to the agents' matrix
        
        history =binornd(1,0.5,1,memory); %Generating first history
        
        storeMinorityGroup = zeros(nbrOfTimeSteps,1);
        storeAttendance = zeros(nbrOfTimeSteps,1); %Resetting
        
        for timeStep = 1:nbrOfTimeSteps
            
            actions = zeros(nbrOfAgents,1);%Resetting
            
            histIndex = bi2de(history)+1;%Histories viewed as a binary
            %number corresponding to a number between 1 and 2^m
            
            
            for i=1:nbrOfAgents  % Which will go
                %Selects the choice perscribed by the at this point most
                %succesful strategy. If there is more than one at the top
                %a coin toss decides the outcome.
                %With a possibility to abstained if the threshold of a
                %stratey is not met.
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
            price(timeStep+1)=price(timeStep)+sum(actions);
            %priceIncr(timeStep)=sum(actions);
            
            
            for i=1:nbrOfAgents %Update score
                
                for j=1:s %for each strategy
                    points(i,j) = points(i,j) - agents(i,histIndex,j)*sum(actions);%-epsilon;
                    %points(i,j) = points(i,j) + agents(i,histIndex,j)*minority;
                end
            end
            
            
            
            storeAttendance(timeStep) = sum(abs(actions)); %Number of agents
            %that doesn't abstain.
            
            %Updates the score. In order to use trick with binary numbers
            %-1 needs to yield a 0. If there is an equal numbers of -1 and
            %1:s due to agents abstaining a 1 or 0 is being coin tossed.
            
            if minority==-1
                history = [history(2:end) 0]; % Update history
                storeMinorityGroup(timeStep)=sum(actions==minority);
                
            elseif minority==0
                storeMinorityGroup(timeStep)=length(find(actions==1));
                history = [history(2:end) binornd(1,0.5,1,1)];
                
            else
                history = [history(2:end) 1]; % Update history
                storeMinorityGroup(timeStep)=sum(actions==minority);
            end
            
            
        end %timestep
        %Store the measure of standard deviation.
        sigma(iRun,memory)=std(storeMinorityGroup);
        
        
    end % End of Runs
    
    
end     %memory

%% Plottar lite goe grejer

figure(1)
clf
subplot(2,1,1)
title('Size of Minority Group')
hold on
plot(storeMinorityGroup)
axis([0 nbrOfTimeSteps 0 nbrOfAgents/2+5])

subplot(2,1,2)
title('Participants')
hold on
plot(storeAttendance)
axis([0 nbrOfTimeSteps 0 nbrOfAgents+5])



dt=1;
for i=1+dt:nbrOfTimeSteps
    logReturn(i)=log(price(i)/price(i-dt));
end
logReturn2=logReturn(1,2:end);
normBrus=normrnd(mean(logReturn2),(std(logReturn2)),1,length(logReturn2));
k=5
figure(2)
clf
subplot(2,1,1)
title('logReturn')
axis([0,length(logReturn2),-k*max(logReturn2),k*max(logReturn2)])
hold on
plot(logReturn2)
subplot(2,1,2)
hold on
title('Normalfördelat brews')
axis([0,length(logReturn2),-k*max(logReturn2),k*max(logReturn2)])
plot(normBrus)

% title('Fördelning av prisökn')
% hist(logReturn2,100)

figure(3)
clf
qqplot(logReturn)







%%
clf

x=linspace(2,12);
figure (1)
axis([0 13 0 15])
hold on
plot(x,0.5*sqrt(nbrOfAgents)*ones(100,1),'--') %Plottar teoretisk
%standardavvikelse för binomialförd.

for k=2:memory
    
    plot(k*ones(nbrOfRuns,1),sigma(:,k),'k.')
    hold on
    
end




