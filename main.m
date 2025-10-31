% Set problem parameters
problem.nvar = 30;
problem.varmin = zeros(1, 30);
problem.varmax = ones(1, 30);
problem.f = @zdt3;

% Set algorithm parameters
params.max_gen = 100;
params.npop = 100;
params.a = 2;
params.a_min = 0;
params.b = 1;
params.elite_ratio = 0.2; %20% Elite Ratio
params.no_improve_threshold = 3;  % Three consecutive generations without improvement triggered the Levy flight.

% RUN IMOWOA
[final_pop, final_obj, pareto_front] = imowoa(problem, params);

% Draw the Pareto frontier
figure;
plot(pareto_front(:,1), pareto_front(:,2), 'ro', 'MarkerSize', 6);
xlabel('f_1');
ylabel('f_2');
title('IMOWOA Pareto frontier');
grid on;