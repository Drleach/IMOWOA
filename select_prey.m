function prey_idx = select_prey(fronts)
    if ~isempty(fronts{1})
        front = fronts{1};
    else
        front = fronts{1};
    end
    prey_idx = front(randi(length(front)));
end