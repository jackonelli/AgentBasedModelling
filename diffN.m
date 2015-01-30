clear all
clf

s = 2;
memory = 2;
slutM=11;
nbrOfTimeSteps = 5000;
nbrOfRuns = 1;

x=linspace(2,12);
figure (1)
axis([0 13 0 15])
hold on
%plot(x,0.5*sqrt(nbrOfAgents)*ones(100,1),'--') %Plottar teoretisk
%standardavvikelse för binomialförd.
    
N=[11 25 101 1001];
z=zeros(length(N),memory-1);
sigma=zeros(length(N),memory-1);


     for a=1:length(N)
        
     nbrOfAgents=N(a);
    
    for memory = 2:slutM

        % Run simulation

        possibleOutcomes = combn([0 1],memory);



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


            sigma=std(storeMinorityGroup);
            varN(a,memory-1)=sigma^2/N(a);
            z(a,memory-1)=2^memory/N(a);

   
    
    end %memory
    
    

     end
%%     
clf
 figure(2)
 marker=['+' '*' 'x' 'o'];
 
 x=linspace(0.0001,10000);
    
    loglog(x,0.25*ones(length(x),1),'--')
    hold on
 for b=2:length(N);
 loglog(z(b,:),varN(b,:),marker(b))
 hold on
 end







