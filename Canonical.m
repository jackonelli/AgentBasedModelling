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
price=ones(nbrOfTimeSteps+1,1)*1000;
priceIncr=zeros(nbrOfTimeSteps,1);

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
          
            histIndex = bi2de(history)+1;
         
            
            for i=1:nbrOfAgents  % Which will go
                               
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
              
                
          %Ren skillnad mellan köp och sälj
            price(timeStep+1)=price(timeStep)+sum(actions);
            %priceIncr(timeStep)=sum(actions);

            
              for i=1:nbrOfAgents %Update score
                  
                  for j=1:s %for each strategy
                      points(i,j) = points(i,j) + agents(i,histIndex,j)*minority;
                      %points(i,j) = points(i,j) + agents(i,histIndex,j)*minority;
                  end
              end
              
            storeMinorityGroup(timeStep)=sum(actions==minority);
            
            storeAttendance(timeStep) = sum(abs(actions)); %Hur många
            %som gör något.
            
            if minority==-1
                history = [history(2:end) 0]; % Update history
            
            elseif minority==0 
                history =binornd(1,0.5,1,memory);
                    
            else
                history = [history(2:end) 1]; % Update history
            end

            
        end %timestep
        
        
        sigma(iRun,memory)=std(storeMinorityGroup);
        
        
    end % End of Runs
    
    
end     %memory

%%

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

%Log-return av pris

dt=20;
for i=1+dt:nbrOfTimeSteps
   logReturn(i)=log(price(i)/price(i-dt));
end

figure(2)
clf
subplot(2,1,1)
title('logReturn')
hold on
plot(logReturn)
subplot(2,1,2)
hold on
title('Fördelning av prisökn')
hist(logReturn,100)

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





