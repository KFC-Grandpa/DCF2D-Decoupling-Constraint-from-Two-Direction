function [ Offspring ] = OperatorDE_pbest_1_main(Population, popsize, Problem, Fitness, p,N)
% The operator of DE


    permutation = randperm(N);
    r0 = permutation;
    [r1,r2,r3] = gnR1R2R3(N, r0);

    array = permutation(1:popsize);
    pop1 = Population(array);

    [~, indBest] = sort(Fitness, 'ascend');
    pNP = max(round(p * N), 2);         % choose at least two best solutions  
    randindex = ceil(rand(1, popsize) * pNP);	% select from [1, 2, 3, ..., pNP]
    randindex = max(1, randindex);              % to avoid the problem that rand = 0 and thus ceil(rand) = 0
    pbest = Population(indBest(randindex));     % randomly choose one of the top 100p% solutions

    Offspring = OperatorDE_pbest_1(Problem,Population(array),pbest,Population(r1(1:popsize)),Population(r2(1:popsize)));
end