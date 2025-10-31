function crowd_dist = crowding_distance1(obj, fronts)
    [n, m] = size(obj);
    crowd_dist = zeros(n, 1);

    for f = 1:length(fronts)
        front = fronts{f};
        front_size = length(front);

        if front_size == 1
            crowd_dist(front) = Inf;
            continue;
        end

        cd = zeros(front_size, 1);
        for m_idx = 1:m
            [~, idx] = sort(obj(front, m_idx));
            sorted_front = front(idx);
            cd(idx(1)) = Inf;
            cd(idx(end)) = Inf;
            f_max = obj(sorted_front(end), m_idx);
            f_min = obj(sorted_front(1), m_idx);
            f_range = f_max - f_min;

            if f_range == 0
                continue;
            end

            for i = 2:front_size-1
                cd(idx(i)) = cd(idx(i)) + (obj(sorted_front(i+1), m_idx) - obj(sorted_front(i-1), m_idx)) / f_range;
            end
        end
        crowd_dist(front) = cd;
    end
end