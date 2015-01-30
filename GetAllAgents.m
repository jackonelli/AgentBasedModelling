function [ agents ] = GetAllAgents(N,memory,s)

agents = [];
m = 2^memory;

    for i=1:N

        %tmpStrat=randi([-1,1]);   
        tmpStrat=binornd(1,0.5,m,s);
                 
        agents = [agents; struct('strategies',tmpStrat)];
    end

end

