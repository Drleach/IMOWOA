function [fronts, rank] = non_dominated_sort(obj)
    [n, m] = size(obj);
    fronts = cell(n, 1);
    rank = zeros(n, 1);
    n_dominated = zeros(n, 1);
    S = cell(n, 1);
    
    for p = 1:n
        for q = 1:n
            if p ~= q
                if all(obj(p,:) <= obj(q,:)) && any(obj(p,:) < obj(q,:))
                    S{p} = [S{p}; q];
                elseif all(obj(q,:) <= obj(p,:)) && any(obj(q,:) < obj(p,:))
                    n_dominated(p) = n_dominated(p) + 1;
                end
            end
        end
        if n_dominated(p) == 0
            fronts{1} = [fronts{1}; p];
            rank(p) = 1;
        end
    end
    
    current_front = 1;
    while ~isempty(fronts{current_front})
        next_front = [];
        for i = 1:length(fronts{current_front})
            p = fronts{current_front}(i);
            for j = 1:length(S{p})
                q = S{p}(j);
                n_dominated(q) = n_dominated(q) - 1;
                if n_dominated(q) == 0
                    next_front = [next_front; q];
                    rank(q) = current_front + 1;
                end
            end
        end
        current_front = current_front + 1;
        fronts{current_front} = next_front;
    end
    
    empty_idx = cellfun(@isempty, fronts);
    fronts(empty_idx) = [];
end