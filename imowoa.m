function [final_pop, final_obj, pareto_front] = imowoa(problem, params)
% Multi-objective Whale Optimization Algorithm Incorporating EMAP Strategy(IMOWOA)
% Input:
%   problem - Problem structure (same as MOWOA)
%   params - Parameter structure, including additional parameters related to EMAP
% Output:
%   final_pop - Final population
%   final_obj - Final population objective values
%   pareto_front - Final Pareto frontier 
% Extract problem parameters
    f = problem.f;
    nvar = problem.nvar;
    varmin = problem.varmin;
    varmax = problem.varmax;
    
% Extract algorithm parameters
    max_gen = params.max_gen;
    npop = params.npop;
    a = params.a;
    a_min = params.a_min;
    b = params.b;
    elite_ratio = params.elite_ratio;  % EMAP elite ratio (e.g. 0.2)
    no_improve_threshold = params.no_improve_threshold;  % The number of times without improvement required to trigger Levi's flight (e.g. 3)
    
% Initialize the population
    pop = initialization(npop, nvar, varmin, varmax);
    
% Calculate the initial target value
    obj = zeros(npop, size(f(pop(1,:)),2));
    for i = 1:npop
        obj(i,:) = f(pop(i,:));
    end
    
% Initialization of EMAP parameters
    elite_num = max(1, floor(npop * elite_ratio));  % Number of elites (at least 1)
    no_improve_count = 0;  % Continuous improvement-free count
    prev_pareto_perf = Inf;  % The previous round's Pareto frontier performance indicators (used for evaluating improvements)
    
   % Main loop
    for gen = 1:max_gen
        % Non-dominated sorting and congestion degree calculation
        [fronts, ~] = non_dominated_sort(obj);
        crowd_dist = crowding_distance1(obj, fronts);
        
        % Update the exploration-development factor 'a' of WOA
        current_a = a - (a - a_min) * (gen / max_gen);
        
        % Generate the offspring population (the core search strategy of WOA)
        pop_offspring = zeros(npop, nvar);
        for i = 1:npop
            % Selecting prey from the non-dominated frontier
            prey_idx = select_prey(fronts);
            prey = pop(prey_idx,:);
            
            % The three search behaviors of WOA
            r1 = rand; r2 = rand;
            A = 2 * current_a * r1 - current_a;
            C = 2 * r2;
            l = rand * 2 - 1;
            p = rand;
            
            if p < 0.5
                if abs(A) < 1
                    D = abs(C * prey - pop(i,:));
                    pop_offspring(i,:) = prey - A * D;
                else
                    rand_idx = randi(npop);
                    X_rand = pop(rand_idx,:);
                    D = abs(C * X_rand - pop(i,:));
                    pop_offspring(i,:) = X_rand - A * D;
                end
            else
                D = abs(prey - pop(i,:));
                pop_offspring(i,:) = D * exp(b * l) * cos(2 * pi * l) + prey;
            end
            
            % Boundary handling
            pop_offspring(i,:) = max(min(pop_offspring(i,:), varmax), varmin);
        end
        
        % Calculate the target value of the offspring
        obj_offspring = zeros(npop, size(obj,2));
        for i = 1:npop
            obj_offspring(i,:) = f(pop_offspring(i,:));
        end
        
% ------------------- Implementation of the EMAP strategy -------------------
% 1. Select the elite population (choose individuals with high crowding degree from the non-dominated frontier)
% Merge the parent and offspring populations for elite selection
        pop_combined = [pop; pop_offspring];
        obj_combined = [obj; obj_offspring];
        [fronts_combined, ~] = non_dominated_sort(obj_combined);
        crowd_dist_combined = crowding_distance1(obj_combined, fronts_combined);
        
        % Select the elites (priority given to individuals with higher congestion levels) from the first non-dominated frontier.
        if ~isempty(fronts_combined{1})
            first_front = fronts_combined{1};
            [~, sorted_idx] = sort(crowd_dist_combined(first_front), 'descend');  % % Descending order of congestion degree
            selected_elite = first_front(sorted_idx(1:min(elite_num, length(first_front))));
        else
      % Extreme case: Selecting elites from all individuals
            [~, sorted_idx] = sort(cellfun(@(x) mean(obj_combined(x,:)), fronts_combined));
            selected_elite = fronts_combined{sorted_idx(1)}(1:min(elite_num, npop));
        end
        elite_positions = pop_combined(selected_elite, :);
        elite_objs = obj_combined(selected_elite, :);
        
% 2. Adaptive selection of perturbation methods (based on the iteration stage and population diversity)
%   Calculate the diversity of the decision space (standard deviation of positions)
        current_diversity = std(pop_combined(:));
        % In the early stage of iteration (the first 50%), exploration is the main focus; in the later stage, development takes precedence.
        if gen < 0.5 * max_gen
            if current_diversity < 1e-5  % Low diversity → Expand search using Levy flight
                perturbation_type = 'levy';
            else  % High diversity → Use uniform perturbation
                perturbation_type = 'uniform';
            end
        else
            % Later stage: If there is no continuous improvement, 
            % a Levy flight will be triggered; otherwise, 
            % a normal perturbation will be used for local development.
            if no_improve_count >= no_improve_threshold
                perturbation_type = 'levy';
            else
                perturbation_type = 'normal';
            end
        end
        
        % 3. Elite individuals generate EMAP offspring
        elso_pop = [];
        elso_obj = [];
        for e = 1:length(selected_elite)
            e_pos = elite_positions(e, :);
            % 生Generate the disturbance term δ 
            switch perturbation_type
                case 'uniform'
                    delta = rand(1, nvar) * 2 - 1;  % U(-1,1)
                case 'normal'
                    delta = randn(1, nvar);  % N(0,1)
                case 'levy'
                    delta = levy_flight(nvar);  % Levi Flight
            end
            
            % Generate EMAP offspring and perform boundary processing
            elso_pos = e_pos + delta .* (varmax - varmin) * 0.1;  % % Scale the disturbance according to the variable range
            elso_pos = max(min(elso_pos, varmax), varmin);
            elso_pop = [elso_pop; elso_pos];
            elso_obj = [elso_obj; f(elso_pos)];
        end
        
        % 4. Merge all individuals (parent generation + offspring generation + EMAP offspring generation)
        pop_all = [pop_combined; elso_pop];
        obj_all = [obj_combined; elso_obj];
        
        % ------------------- Environmental selection (maintaining population size)） -------------------
        [fronts_all, ~] = non_dominated_sort(obj_all);
        crowd_dist_all = crowding_distance1(obj_all, fronts_all);
        [pop, obj] = environmental_selection(pop_all, obj_all, fronts_all, crowd_dist_all, npop);
        
% 5. Update the count of no-improvement cases (based on Pareto frontier performance)
%   Use the average target value of the first frontier to evaluate the improvements (a simplified approach in the multi-objective scenario)
        current_fronts = non_dominated_sort(obj);
        current_pareto = obj(current_fronts{1}, :);
        current_perf = mean(current_pareto(:));  % Average target value (the smaller, the better)
        
        if current_perf < prev_pareto_perf - 1e-6  % There are improvements (with a small threshold to avoid floating-point errors)
            no_improve_count = 0;
            prev_pareto_perf = current_perf;
        else
            no_improve_count = no_improve_count + 1;
        end
        
        % Display iteration information
        fprintf('Generation %d/%d, No-improve count: %d\n', gen, max_gen, no_improve_count);
    end
    
    % Extract the final Pareto frontier
    [fronts, ~] = non_dominated_sort(obj);
    pareto_front = obj(fronts{1},:);
    final_pop = pop;
    final_obj = obj;
end








function [pop, obj] = environmental_selection(pop_combined, obj_combined, ...
        fronts_combined, crowd_dist_combined, npop)
    pop = [];
    obj = [];
    remaining = npop;
    
    for f = 1:length(fronts_combined)
        front = fronts_combined{f};
        front_size = length(front);
        
        if remaining >= front_size
            pop = [pop; pop_combined(front,:)];
            obj = [obj; obj_combined(front,:)];
            remaining = remaining - front_size;
        else
            [~, idx] = sort(crowd_dist_combined(front), 'descend');
            selected = front(idx(1:remaining));
            pop = [pop; pop_combined(selected,:)];
            obj = [obj; obj_combined(selected,:)];
            remaining = 0;
            break;
        end
    end
end



