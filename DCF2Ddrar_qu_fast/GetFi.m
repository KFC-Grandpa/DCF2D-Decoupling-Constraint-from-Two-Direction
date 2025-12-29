function [NewArchive,Fitness,flag] = GetFi(Pop1,Archive,Offspring,N)
% GetFi: Update Infeasible Archive (LI) and return update status
% Pop1: Main Population (MT)
% Archive: Current Infeasible Archive (Old LI)
% Offspring: Newly generated offspring
% N: Archive Size limit

    % 1. Combine Old Archive and Offspring to form candidates
    Pop2 = [Archive, Offspring];
    
    Pop1Obj = Pop1.objs;
    Pop2Obj = Pop2.objs;
    
    N1 = size(Pop1Obj,1);
    N2 = size(Pop2Obj,1);

    %% 2. Filter 1: Retain solutions that dominate MT but are not dominated by MT
    % Dominate(i): How many solutions in MT are dominated by Candidate i
    % save(i): False if Candidate i is dominated by any solution in MT
    Dominate = zeros(1,N2);
    save = true(1,N2);
    
    for i = 1 : N2
        for j = 1 : N1
            % k = 1: Candidate i dominates MT j
            % k = -1: MT j dominates Candidate i
            k = any(Pop2Obj(i,:)<Pop1Obj(j,:)) - any(Pop2Obj(i,:)>Pop1Obj(j,:));
            
            if k == -1 
                save(i) = false; % Candidate dominated by MT -> Discard
            elseif k == 1
                Dominate(i) = Dominate(i) + 1; % Candidate dominates MT -> Good
            end
        end
    end
    
    % Only keep solutions that dominate at least one MT solution AND are not dominated
    save(Dominate==0) = 0; 
    Pop2 = Pop2(save);
    
    %% 3. Filter 2: Environmental Selection
    % Priority: Infeasible > Feasible
    % For Infeasible: Use Negative Objectives (Closer to CPF boundary)
    if length(Pop2) > N
        Pop2Con = Pop2.cons;
        Pop2Con = max(0,Pop2Con);
        CV = sum(Pop2Con,2);
        
        % Set CV of feasible solutions to huge value to put them in the last front
        CV(CV==0) = 999999;
        % Restore CV for infeasible ones (strictly > 0)
        CV(CV<999999 & CV>0) = 0; % Treat all real infeasible as same "layer" for NDSort?
        % Note: Your original code had CV(CV<999999)=0. 
        % This effectively makes NDSort rank purely based on -objs for all infeasible solutions.
        % This is correct for the logic "Find closest to CPF from infeasible side".
        CV(CV<999999) = 0; 
        
        % Sort based on Negative Objectives (Maximize Objectives -> Closer to CPF)
        FrontNo = NDSort(-Pop2.objs,CV,inf);
        
        % Calculate Crowding Distance (Standard SPEA2/NSGAII logic)
        Pop2Obj = Pop2.objs;
        Distance = pdist2(Pop2Obj,Pop2Obj);
        Distance(logical(eye(length(Distance)))) = inf;
        Distance = sort(Distance,2);
        D = 1./(Distance(:,floor(sqrt(length(Pop2)))+2));
        
        Fit = FrontNo + D';
        [Fit,Fitnessind] = sort(Fit);
        
        NewArchive = Pop2(Fitnessind(1:N));
        Fitness = Fit(1:N);
    else
        NewArchive = Pop2;
        Fitness = []; % Calculate if needed or return empty
    end
    
    %% 4. Calculate Flag (Has the Archive Changed?)
    if isempty(Archive)
        flag = ~isempty(NewArchive);
    elseif length(Archive) ~= length(NewArchive)
        flag = true;
    else
        % Use sortrows to ignore order differences
        % If objective values are identical, we consider the archive unchanged
        flag = ~isequal(sortrows(Archive.objs), sortrows(NewArchive.objs));
    end
end