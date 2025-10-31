function [final_pop, final_obj, pareto_front] = mowoa(problem, params)

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
    
    % Initialize the population
    pop = initialization(npop, nvar, varmin, varmax);
    
    % Calculate the target value of the initial population
    obj = zeros(npop, size(f(pop(1,:)),2));  % Automatically adjust to the target quantity
    for i = 1:npop
        obj(i,:) = f(pop(i,:));
    end
    
    % Main loop
    for gen = 1:max_gen
        % Calculate the non-dominated sorting and crowding degree of the current population
        [fronts, ~] = non_dominated_sort(obj);
        crowd_dist = crowding_distance(obj, fronts);
        
        % Update exploration - development factor a (linearly decreasing)
        current_a = a - (a - a_min) * (gen / max_gen);
        
        % Generate the offspring population
        pop_offspring = zeros(npop, nvar);
        for i = 1:npop
            % Randomly select an individual from the non-dominated frontier as the "prey" (reference solution)
            prey_idx = select_prey(fronts);
            prey = pop(prey_idx,:);
            
      
            r1 = rand;  % [0,1]
            r2 = rand;  % [0,1]
            A = 2 * current_a * r1 - current_a;  
            C = 2 * r2;                          
            l = rand * 2 - 1;                     % [-1,1]
            p = rand;                             % [0,1]
            
            % Three search strategies of the whale optimization algorithm
            if p < 0.5
                if abs(A) < 1
                    % Surround the prey
                    D = abs(C * prey - pop(i,:));
                    pop_offspring(i,:) = prey - A * D;
                else
                    % 2. Random search (selecting random individuals from the population)
                    rand_idx = randi(npop);
                    X_rand = pop(rand_idx,:);
                    D = abs(C * X_rand - pop(i,:));
                    pop_offspring(i,:) = X_rand - A * D;
                end
            else
                % 3. Bubble Network Attack (Spiral Update)
                D = abs(prey - pop(i,:));
                pop_offspring(i,:) = D * exp(b * l) * cos(2 * pi * l) + prey;
            end
            
            % Boundary handling (ensuring that the variables remain within the feasible region)
            pop_offspring(i,:) = max(min(pop_offspring(i,:), varmax), varmin);
        end
        
        % Calculate the target value of the offspring population
        obj_offspring = zeros(npop, size(obj,2));
        for i = 1:npop
            obj_offspring(i,:) = f(pop_offspring(i,:));
        end
        
        % Merge the parental and offspring populations
        pop_combined = [pop; pop_offspring];
        obj_combined = [obj; obj_offspring];
        
        % Perform non-dominated sorting and crowding degree calculation for the merged population
        [fronts_combined, ~] = non_dominated_sort(obj_combined);
        crowd_dist_combined = crowding_distance(obj_combined, fronts_combined);
        
        % Environmental selection (retaining the superior individuals and maintaining the population size)
        [pop, obj] = environmental_selection(pop_combined, obj_combined, ...
            fronts_combined, crowd_dist_combined, npop);
        
        % Display current iteration information
        fprintf('Generation %d/%d completed\n', gen, max_gen);
    end
    
    % Extract the final Pareto frontier (the first non-dominated layer)
    [fronts, ~] = non_dominated_sort(obj);
    pareto_front = obj(fronts{1},:);
    final_pop = pop;
    final_obj = obj;
end







