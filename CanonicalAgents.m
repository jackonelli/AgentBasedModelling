function [agents] = CanonicalAgents(N,P,nbrOfStrategies)
agents=ones(N,P,nbrOfStrategies)-floor(rand(N,P,nbrOfStrategies)*2)*2;
end
