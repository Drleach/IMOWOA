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