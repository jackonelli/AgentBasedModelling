function [agents] = CanonicalAgents(N,memory,nbrOfStrategies)
P = 2^memory;
agents=ones(N,P,nbrOfStrategies)-floor(rand(N,P,nbrOfStrategies)*2)*2;
end
