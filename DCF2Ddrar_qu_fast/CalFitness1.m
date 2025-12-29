function Fitness = CalFitness1(Pop2Obj,Pop1Obj,Pop2Con)
% Calculate the fitness of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    N1 = size(Pop1Obj,1);
    N2 = size(Pop2Obj,1);
    CV = sum(max(0,Pop2Con),2);
    %% Detect the dominance relation between each two solutions
    Dominate = zeros(1,N2);
    for i = 1 : N2
        for j = 1 : N1
            k = any(Pop2Obj(i,:)<Pop1Obj(j,:)) - any(Pop2Obj(i,:)>Pop1Obj(j,:));
            %1支配2k=-1
            if k == -1 
                Dominate(i) = Dominate(i) + 1;  
            end
        end
    end
    Betterindex = find(Dominate==0);
    BetterObj = Pop2Obj(Betterindex,:);
    BetterCon = Pop2Con(Betterindex,:);

    
    if size(BetterObj,1) > N1
        Badindex = find(Dominate>0);
        BadObj = Pop2Obj(Badindex,:);
        BadCon = Pop2Obj(Badindex,:);
        FitBad = CalFitness(0,BadObj,BadCon);
        FitBet = CalFitness(0,-BetterObj,-BetterCon);

        Fitness = zeros(N2,1);
        Fitness(Badindex) = FitBad;
        Fitness(Betterindex) = FitBet;

    else
        Distance = pdist2(Pop2Obj,Pop2Obj);
        Distance(logical(eye(length(Distance)))) = inf;
        Distance = sort(Distance,2);
        D = 1./(Distance(:,floor(sqrt(N2)))+2);

        Fitness = Dominate + D';
    end
    
    % Do2 = Pop2Obj(Dominate);
    % 
    % if size(Do2) <= N1
    %     Fitness = true(1,N2);
    %     Fitness(Dominate==1) = 0;
    % else
    %% Calculate D(i)
    

    
end