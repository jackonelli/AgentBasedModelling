function [ agents ] = CanonicalAgents(N,memory,s)

m = 2^memory;
agents = zeros(N,m,s);
    for i=1:N

        agents(i,:,:) =ones(m,s) - floor(rand(m,s)*2)*2;
        
    end
    
end
