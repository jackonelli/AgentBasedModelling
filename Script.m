clear all
clf

nbrOfAgents = 101;
s = 2;
memory = 2;
slutM=11;
nbrOfTimeSteps = 500;
nbrOfRuns = 10;
sigma=ones(nbrOfRuns,slutM);

x=linspace(2,12);
figure (1)
axis([0 13 0 15])
hold on
plot(x,0.5*sqrt(nbrOfAgents)*ones(100,1),'--') %Plottar teoretisk
%standardavvikelse för binomialförd.



for memory = 2:slutM
    
    % Run simulation
    
    possibleOutcomes = combn([0 1],memory);
    
    for iRun = 1:nbrOfRuns
        
        agents = GetAllAgents(nbrOfAgents,memory,s);
        points = zeros(nbrOfAgents, s);
        history0 = [zeros(1,floor(memory/2)) ones(1,memory-floor(memory/2))];
        history = history0(randperm(numel(history0)));
        storeMinorityGroup = zeros(nbrOfTimeSteps,1);
        
        for timeStep = 1:nbrOfTimeSteps
            attendance = zeros(nbrOfAgents,1); %Varför nollställa här?
            %storeIndex=zeros(nbrOfAgents,s);
            tmpStrategies = [];
            histIndex=-1;
            
            for j=1:size(possibleOutcomes,1)
                if isequal(possibleOutcomes(j,:),history)
                    histIndex = j;
                    break
                end
            end
            
            for i=1:nbrOfAgents  % Which will go
                
                if points(i,1) == points(i,2)
                    index = floor(rand*2)+1;
                else
                    [~,index] = max(points(i,:));
                end
                %index = index-s; % Best strategy
                tmpStrategies = agents(i).strategies;
                tmpStrategy = tmpStrategies(:,index);
                attendance(i) = tmpStrategy(histIndex);
                
            end
            
            minority = sum(attendance) < nbrOfAgents/2;
            
            for i=1:nbrOfAgents %Update score
                tmpStrategies = agents(i).strategies;
                for j=1:s %for each strategy
                    points(i,j) = points(i,j) + (tmpStrategies(histIndex,j) == minority);
                end
                               
            end
            
            storeMinorityGroup(timeStep) = sum(attendance);%sum(attendance == minority);
            %För att få det att stämma överens måste man ha bara
            %attendance och inte storleken på minoritetsgruppen.
            
            history = [history(2:end) minority]; % Update history
            
        end %timestep
        
        
        sigma(iRun,memory)=std(storeMinorityGroup);
        
        
    end % End of Runs
    
    
end %memory


for k=2:memory
    
    plot(k*ones(nbrOfRuns,1),sigma(:,k),'k.')
    hold on
    
end









