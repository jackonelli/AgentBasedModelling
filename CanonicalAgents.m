function [agents] = CanonicalAgents(N,memory,nbrOfStrategies)

m = 2^memory;
agents = zeros(N,m,nbrOfStrategies);
%points = (epsilon+1)*ones(nbrOfAgents, nbrOfStrategies);

    for i=1:N

        agents(i,:,:) =ones(m,nbrOfStrategies) - floor(rand(m,nbrOfStrategies)*2)*2;
        
    end
    
end
