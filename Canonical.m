clear all

nbrOfAgents = 101;
s = 2;
slutM=6;
nbrOfTimeSteps = 1000;
nbrOfRuns = 1;
sigma=ones(nbrOfRuns,slutM);
kvot=0.46;
counter=0;
lambda=1.5;
price=ones(nbrOfTimeSteps+1,1);


for memory = 2:slutM
    
    % Run simulation
    
    %possibleOutcomes = combn([0 1],memory);
    
    for iRun = 1:nbrOfRuns
      
        agents = CanonicalAgents(nbrOfAgents,memory,s);
       
        points = ones(nbrOfAgents, s);
        
        history =binornd(1,0.5,1,memory); %Genererar
        %första historia.
   
        storeMinorityGroup = zeros(nbrOfTimeSteps,1);
        storeAttendance = zeros(nbrOfTimeSteps,1);
        
        for timeStep = 1:nbrOfTimeSteps
            actions = zeros(nbrOfAgents,1); 
            %storeIndex=zeros(nbrOfAgents,s);
            tmpStrategies = [];
            histIndex=-1;
            
            %for j=1:size(possibleOutcomes,1)
             %   if isequal(possibleOutcomes(j,:),history)
                    histIndex = bi2de(history)+1;
              %      break
               % end
            %end
            
            for i=1:nbrOfAgents  % Which will go
                %Per ska ändra till generella stratIndex
                
                if max(points(i,:))>kvot*timeStep
                    
                  
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
            
            if sum(actions)==0
            counter=counter+1;
            end
            
            
            
            minority=-sign(sum(actions));
            price(timeStep+1)=price(timeStep)*exp(sum(actions)/lambda);
            
            for i=1:nbrOfAgents %Update score
                
                for j=1:s %for each strategy
                      points(i,j) = points(i,j) + (agents(i,histIndex,j)== minority);
                      %points(i,j) = points(i,j) + agents(i,histIndex,j)*minority;
                end                           
                
                               
            end
            
            storeMinorityGroup(timeStep)=sum(actions==minority);
            
            storeAttendance(timeStep) = sum(abs(actions)); %Hur många
            %som gör något.
            
            if minority==-1
            history = [history(2:end) 0]; % Update history
            
            else
            history = [history(2:end) 1]; % Update history
            end

            
        end %timestep
        
        
        sigma(iRun,memory)=std(storeMinorityGroup);
        
        
    end % End of Runs
    
    
end     %memory

%%
clf
figure(2)

subplot(3,1,1)
title('Size of Minority Group')
hold on
plot(storeMinorityGroup)
axis([0 nbrOfTimeSteps 0 nbrOfAgents/2+5])

subplot(3,1,2)
title('Participants')
hold on
plot(storeAttendance)
axis([0 nbrOfTimeSteps 0 nbrOfAgents+5])

figure(3)
clf
title('log(Price)')
%hold on
semilogy(price)




frequency=counter/nbrOfTimeSteps;




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





