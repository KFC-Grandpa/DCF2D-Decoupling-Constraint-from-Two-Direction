classdef DCF2Ddrar_qu_fast < ALGORITHM
% <2025> <multi> <real/integer/label/binary/permutation> <constrained>
% DCF2D: Decoupling Constraint from Two Directions

    methods
        function main(Algorithm,Problem)
            %% ================= 1. Initialization =================
            gen = 2;
            
            % --- Main Population (MT) ---
            Population1 = Problem.Initialization(); 
            Fitness1    = CalFitness(0,Population1.objs,Population1.cons);
            Zmin1       = min(Population1.objs,[],1);
            Lastinfes=[];
            N_AT = floor(Problem.N/4);

            % --- Infeasible Archive (LI) ---
            [Lastinfes, ~, ~] = GetFi(Population1, Lastinfes, Population1, Problem.N);
            
            % --- Auxiliary Tasks (ATs) ---
            % Pops{1...totalcon}: Constraint-specific ATs
            % Pops{totalcon+1}:   UPF Search AT (AT_0)
            transtage = 0; % 0: Stage 1, 1: Stage 2
            totalcon = size(Population1.cons,2);

            for i = 1 : totalcon+1
                if i == totalcon+1
                    Pops{i}.N = Problem.N;
                else
                    Pops{i}.N = N_AT;
                end
                Pops{i}.pop = Population1(1:Pops{i}.N);
                Pops{i}.fit = CalFitness(0, Pops{i}.pop.objs);
                Pops{i}.dir = 1; % 1: Positive (SCPF), 0: Negative (RCPF)
                Pops{i}.Z = min(Pops{i}.pop.objs,[],1);
                Pops{i}.act = 0;
                Pops{i}.color = [rand,rand,rand]; % For visualization
                Pops{i}.protect_gen=0;
            end

            % Activate AT_0 initially
            Pops{totalcon+1}.act = 1; 
            
            % --- State Variables ---
            cnt = Problem.maxFE; % Counter for periodic detection
            dw = 0; % Draw flag (debug)
            Candidate = [];
            last_Candidate = [];
            
            % --- [NEW] Sliding Window Trend Analysis (Replaces is_stable) ---
            LastCenter_AT0 = inf(1, Problem.M);
            WindowSize = 5; % Observe trend over 5 generations
            MoveHistory_AT0 = ones(1, WindowSize) * 1.0; % Init with high variation

            %% ================= 2. Main Loop =================
            while Algorithm.NotTerminated(Population1)
                
                FE(gen) = Problem.FE;
                Offspring = [];

                %% -------- A. Resource Allocation (Scheme 1) --------
                % Strategy: MT takes N/2, ATs share the rest.
                
                % 1. Count active Auxiliary Tasks
                active_AT_idx = find(cellfun(@(x) x.act, Pops));
                num_active_AT = length(active_AT_idx);

                % 2. Allocate Quotas
                Total_Quota = Problem.N;
                
                % MT gets fixed 50%
                Quota_MT = floor(0.5 * Total_Quota); 
                
                % ATs share the remaining 50%
                if num_active_AT > 0
                    Remaining_Quota = Total_Quota - Quota_MT;
                    Quota_Per_AT = floor(Remaining_Quota / num_active_AT);
                    
                    % [Minimum Guarantee]: DE needs at least 4 individuals (1 target + 3 parents)
                    % We set 5 to be safe. This might slightly exceed Total_Quota but ensures validity.
                    Quota_Per_AT = max(5, Quota_Per_AT); 
                else
                    Quota_Per_AT = 0;
                end

                %% -------- B. Reproduction (Main Task) --------
                MatingPool = [Population1(randsample(Problem.N, Problem.N))];
                [Mate1, Mate2, Mate3] = Neighbor_Pairing_Strategy(MatingPool, Zmin1);
                
                % Split quota between rand and pbest operators
                num_rand_mt = floor(Quota_MT / 2);
                num_pbest_mt = Quota_MT - num_rand_mt;
                
                Off1 = [];
                if num_rand_mt > 0
                    Off1 = [Off1, OperatorDE_rand_1(Problem, Mate1(1:num_rand_mt), Mate2(1:num_rand_mt), Mate3(1:num_rand_mt))];
                end
                if num_pbest_mt > 0
                    % Assumes OperatorDE_pbest_1_main generates 'num_pbest_mt' offspring
                    temp_off = OperatorDE_pbest_1_main(Population1, num_pbest_mt, Problem, Fitness1, 0.1, Problem.N);
                    Off1 = [Off1, temp_off];
                end
                
                if ~isempty(Off1)
                    Zmin1 = min([Zmin1; Off1.objs], [], 1);
                    Offspring = [Offspring, Off1];
                end

                %% -------- C. Reproduction (Auxiliary Tasks) --------
                for i = 1 : totalcon+1 
                    if Pops{i}.act
                        MatingPool = [Pops{i}.pop(randsample(Pops{i}.N, Pops{i}.N))];
                        [Mate1, Mate2, Mate3] = Neighbor_Pairing_Strategy(MatingPool, Pops{i}.Z);
                        
                        num_rand_at = floor(Quota_Per_AT / 2);
                        num_pbest_at = Quota_Per_AT - num_rand_at;
                        
                        Off = [];
                        if num_rand_at > 0
                            Off = [Off, OperatorDE(Problem, Mate1(1:num_rand_at), Mate2(1:num_rand_at), Mate3(1:num_rand_at))];
                        end
                        if num_pbest_at > 0
                            temp_off = OperatorDE_pbest_1_main(Pops{i}.pop, num_pbest_at, Problem, Pops{i}.fit, 0.1, Pops{i}.N);
                            Off = [Off, temp_off];
                        end
                        
                        if ~isempty(Off)
                            Pops{i}.Z = min([Pops{i}.Z; Off.objs], [], 1);
                            Offspring = [Offspring, Off];
                        end
                    end
                end

                %% -------- D. Environmental Selection --------
                % 1. Update MT (Main Population)
                [Population1, Fitness1] = EnvironmentalSelection([Population1, Offspring], Problem.N, true);
                
                % 2. Update LI (Infeasible Archive)
                [Lastinfes, ~, is_LI_updated] = GetFi(Population1, Lastinfes, Offspring, Problem.N);

                if is_LI_updated && transtage == 1 % 只在 Stage 2 生效
                    % 找出当前 LI 中所有解违反的约束索引
                    % max(0, cons) > 0 表示违反，sum(...) > 0 表示该约束至少被一个解违反
                    Violated_Cons_Idx = find(sum(max(0, Lastinfes.cons), 1) > 0);
                    
                    for k = 1 : length(Violated_Cons_Idx)
                        con_id = Violated_Cons_Idx(k);
                        % 如果该约束对应的任务处于“关闭”状态，且之前被误判了
                        if Pops{con_id}.act == 0
                            % 复活它！
                            Pops{con_id}.act = 1;
                            % 并且重置为正向搜索 (或者保持原有策略，通常重置为正向比较稳妥)
                            Pops{con_id}.dir = 1; 

                            Pops{con_id}.protect_gen=5;
                            
                            % 可选：从 LI 中借几个解给它作为种子，帮助它快速启动
                            % (这一步不是必须的，因为下一代它会利用共享的 Offspring)
                        end
                    end
                end
                
                % Visualization (Optional)
                % if dw
                %     Draw(Lastinfes.objs,'*','Markeredgecolor',[.2,.2,.2],'Markerfacecolor',[.1,.4,.1]);
                % end

                % 3. Update ATs (Auxiliary Tasks)
                for i = 1 : totalcon+1
                    if Pops{i}.act
                        if i == totalcon+1 % AT_0 (UPF Task)
                            if ~transtage
                                % Stage 1: Pure UPF search
                                [Pops{i}.pop, Pops{i}.fit] = EnvironmentalSelection([Pops{i}.pop, Offspring], Problem.N, false);
                            else
                                % Stage 2: Explore region between UPF and MT (Negative Direction)
                                [Pops{i}.pop, Pops{i}.fit] = ESConNeg(Population1, [Pops{i}.pop, Offspring], Problem.N, 0);
                            end
                            if dw; Draw(Pops{i}.pop.objs,'o','Markeredgecolor',Pops{i}.color,'Markerfacecolor',Pops{i}.color); end
                        elseif Pops{i}.dir % Positive Direction (SCPF)
                            [Pops{i}.pop, Pops{i}.fit] = ESCon([Pops{i}.pop, Offspring], Pops{i}.N, i);
                            if dw; Draw(Pops{i}.pop.objs,'+','Markeredgecolor',Pops{i}.color,'Markerfacecolor',Pops{i}.color); end
                        else % Negative Direction (RCPF)
                            [Pops{i}.pop, Pops{i}.fit] = ESConNeg(Population1, [Pops{i}.pop, Offspring], Pops{i}.N, i);
                            if dw; Draw(Pops{i}.pop.objs,'sk','Markeredgecolor',Pops{i}.color,'Markerfacecolor',Pops{i}.color); end
                        end
                        if Pops{i}.protect_gen >=0; Pops{i}.protect_gen = Pops{i}.protect_gen - 1; end
                    end
                end

                %% -------- E. Stage Control Logic --------
                
                % === Stage 1 Check: Sliding Window Trajectory Analysis ===
                if transtage == 0 
                    Pop_AT0 = Pops{totalcon+1}.pop;
                    CurCenter_AT0 = mean(Pop_AT0.objs, 1);
                    
                    % 1. Calculate Relative Movement
                    if isinf(LastCenter_AT0(1))
                        RelativeMove = 1;
                    else
                        MoveDist = norm(CurCenter_AT0 - LastCenter_AT0);
                        CenterNorm = norm(CurCenter_AT0) + 1e-10; % Avoid div by zero
                        RelativeMove = MoveDist / CenterNorm;
                    end
                    
                    % 2. Update Sliding Window History
                    MoveHistory_AT0(1:end-1) = MoveHistory_AT0(2:end); % Shift left
                    MoveHistory_AT0(end) = RelativeMove;               % Append new
                    
                    LastCenter_AT0 = CurCenter_AT0;
                    
                    % 3. Check Convergence Trend
                    AvgMove = mean(MoveHistory_AT0);
                    % Threshold 1e-4 implies <0.01% change on average over 5 generations
                    if AvgMove < 1e-4 && gen > WindowSize 
                        transtage = 1; % Switch to Stage 2
                        
                        % 【新增】Stage 2 初始化：立即根据当前 LI 激活任务
                        Pops{totalcon+1}.act = 0; % 关闭 AT_0
                        
                        % 扫描 LI
                        Violated_Cons_Idx = find(sum(max(0, Lastinfes.cons), 1) > 0);
                        for k = 1 : length(Violated_Cons_Idx)
                            con_id = Violated_Cons_Idx(k);
                            Pops{con_id}.act = 1;
                            Pops{con_id}.dir = 1;
                            Pops{con_id}.protect_gen = 5;
                        end
                    end
                    
                else

                    % Dynamic Direction Adjustment
                    for i = 1 : totalcon
                        if Pops{i}.act
                            if Pops{i}.dir 
                                % Positive -> Negative if blocked
                                if min(sum(max(0, Pops{i}.pop.cons), 2)) > 0
                                    Pops{i}.dir = 0;
                                end
                            elseif ~Pops{i}.dir
                                % Negative -> Deactivate if dominated by MT
                                temp = [Pops{i}.pop, Population1];
                                FrontNo = NDSort(temp.objs, inf);
                                % If no solution from AT (indices 1:N) is in the first front
                                if sum(FrontNo(1:Pops{i}.N) == 1) == 0 && Pops{i}.protect_gen == 0
                                     Pops{i}.act = 0;
                                end
                            end
                        end
                    end
                end
                
                % === Stage 3 Check: Final Convergence ===
                if Problem.FE > Problem.maxFE * 0.9
                    for i = 1 : totalcon+1
                        Pops{i}.act = 0;
                    end
                end

                gen = gen + 1;
            end
        end
    end
end