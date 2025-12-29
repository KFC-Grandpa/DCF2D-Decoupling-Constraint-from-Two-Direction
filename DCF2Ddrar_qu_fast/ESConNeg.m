function [return_pop,return_Fitness] = ESConNeg(Pop1,Population,N,processcon)

if processcon
    input_cons = Population.cons;
    if processcon
        input_cons = -input_cons(:,processcon);
    else
        input_cons = -input_cons;
    end
    input_cons(input_cons<0) = 0;
    input_cons = sum(input_cons,2);

    findex = find(input_cons<=0);
    ifindex = find(input_cons>0);

    fPopulation = Population(findex);
    ifPopulation = Population(ifindex);
else
    fPopulation = Population;
end

    if isempty(fPopulation)
        ifFitness = CalFitness(ifPopulation.objs,ifPopulation.cons);
        Next2 = ifFitness < 1;
        if sum(Next2) <= N
            [~,Rank] = sort(ifFitness);
            Next2(Rank(1:N )) = true;
        elseif sum(Next2) > N
            %% Modification Start: Optimize truncation using existing Fitness density
            % Old logic: Del = Truncation(ifPopulation(Next2).objs, sum(Next2)-N ); ...
            
            % Get indices of candidates (Next2 is logical true)
            CandidateIdx = find(Next2);
            % Sort these candidates based on their Fitness (Ascending: Lower fitness = Better sparsity)
            [~, SortIdx] = sort(ifFitness(CandidateIdx));
            % Reset Next2 and only keep the top N
            Next2(:) = false;
            Next2(CandidateIdx(SortIdx(1:N))) = true;
            %% Modification End
        end

        ifPopulation = ifPopulation(Next2);
        ifFitness    = ifFitness(Next2);
        % Sort the population
        [ifFitness,rank] = sort(ifFitness);
        ifPopulation = ifPopulation(rank);

        fPopulation = [];
        fFitness = [];

    elseif length(fPopulation) <= N
        cons = fPopulation.cons;
        cons(cons<0)=0;
        cons = sum(cons,2);
        fFitness = CalFitness(0,[fPopulation.objs,cons]);
        Next = fFitness < 1;

        [~,Rank] = sort(fFitness);
        Next(Rank(1:length(fPopulation) )) = true;

        fPopulation = fPopulation(Next);
        fFitness    = fFitness(Next);
        % Sort the population
        [fFitness,rank] = sort(fFitness);
        fPopulation = fPopulation(rank);
        %
        ifFitness = CalFitness(ifPopulation.objs,ifPopulation.cons); % ,
        Next2 = ifFitness < 1;
        
        NumNeeded = N - length(fPopulation); % Define needed count
        
        if sum(Next2) <= NumNeeded
            [~,Rank] = sort(ifFitness);
            Next2(Rank(1:NumNeeded )) = true;
        elseif sum(Next2) > NumNeeded
            %% Modification Start: Optimize truncation using existing Fitness density
            % Old logic: Del = Truncation(ifPopulation(Next2).objs, sum(Next2)-(N - length(fPopulation)) ); ...
            
            CandidateIdx = find(Next2);
            [~, SortIdx] = sort(ifFitness(CandidateIdx));
            Next2(:) = false;
            Next2(CandidateIdx(SortIdx(1:NumNeeded))) = true;
            %% Modification End
        end

        ifPopulation = ifPopulation(Next2);
        ifFitness    = ifFitness(Next2) + max(fFitness);
        % Sort the population
        [ifFitness,rank] = sort(ifFitness);
        ifPopulation = ifPopulation(rank);

    elseif length(fPopulation) > N

        Pop1Obj = Pop1.objs;
        Pop2Obj = fPopulation.objs;


        N1 = size(Pop1Obj,1);
        N2 = size(Pop2Obj,1);
        discard = false(1,N2);
        for i = 1 : N2
            for j = 1 : N1
                k = any(Pop2Obj(i,:)<Pop1Obj(j,:)) - any(Pop2Obj(i,:)>Pop1Obj(j,:));
                %1支配2k=-1
                if k == -1 
                    discard(i) = true;
                end
            end
        end
        cons = fPopulation.cons;
        cons(cons<0)=0;
        cons = sum(cons,2);
        fFitness = CalFitness(0,[fPopulation.objs,cons]);
        maxfFitness = max(fFitness);
        fFitness(discard) = fFitness(discard) + maxfFitness;
        Next = fFitness < 1;
        if sum(Next) <= N
            [~,Rank] = sort(fFitness);
            Next(Rank(1:N )) = true;
            fPopulation = fPopulation(Next);
            fFitness   = fFitness(Next);
        elseif sum(Next) > N
            if processcon
                NfPopulation = fPopulation(Next);
                NfFitness = CalFitness(0,-[NfPopulation.objs]);
                NextN = NfFitness < 1;
                if sum(NextN) < N
                    [~,Rank] = sort(NfFitness);
                    NextN(Rank(1:N)) = true;
                elseif sum(NextN) > N
                    %% Modification Start: Optimize truncation using NfFitness
                    % Old logic: Del = Truncation(NfPopulation(NextN).objs,sum(NextN)-N); ...
                    
                    CandidateIdx = find(NextN);
                    [~, SortIdx] = sort(NfFitness(CandidateIdx));
                    NextN(:) = false;
                    NextN(CandidateIdx(SortIdx(1:N))) = true;
                    %% Modification End
                end
                fPopulation = NfPopulation(NextN);
                fFitness    = NfFitness(NextN);
            else
                %% Modification Start: Optimize truncation using fFitness
                % Old logic: Del = Truncation(fPopulation(Next).objs, sum(Next)-N ); ...
                
                CandidateIdx = find(Next);
                [~, SortIdx] = sort(fFitness(CandidateIdx));
                Next(:) = false;
                Next(CandidateIdx(SortIdx(1:N))) = true;
                
                fPopulation = fPopulation(Next);
                fFitness   = fFitness(Next);
                %% Modification End
            end
        end

        
        % Sort the population
        [fFitness,rank] = sort(fFitness);
        fPopulation = fPopulation(rank);

        ifPopulation = [];
        ifFitness = [];
    end

    return_pop = [fPopulation,ifPopulation];
    return_Fitness = [fFitness,ifFitness];
end

% Function Truncation is removed/unused as per requirement