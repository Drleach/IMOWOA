% MOWOA unique function: Select "prey" (reference solutions) from the non-dominated frontier
function prey_idx = select_prey(fronts)
    %优先从第一非支配层进行选择，以提高收敛速度。
    if ~isempty(fronts{1})
        front = fronts{1};
    else
        front = fronts{1};  % % If the first layer is empty (the initial generation may be the case), select from any layer
    end
    % Randomly select a non-dominated individual as the prey
    prey_idx = front(randi(length(front)));
end