function [Population, Fitness] = ESCon(Population, N, con_idx)
% Positive Direction: Search for SCPF

    PopObj = Population.objs;
    PopCons = Population.cons;
    
    if con_idx > 0
        current_cv = max(0, PopCons(:, con_idx));
    else

        current_cv = zeros(size(PopObj, 1), 1);
    end
    

    Fitness = CalFitness(0, PopObj, current_cv);
    
    [Fitness, Rank] = sort(Fitness);
    
    if length(Population) > N
        Population = Population(Rank(1:N));
        Fitness = Fitness(1:N);
    else
        Population = Population(Rank);
    end
end