% % Set problem parameters
problem.nvar = 30;                
problem.varmin = zeros(1, 30);    
problem.varmax = ones(1, 30);    
problem.f = @zdt1;                  

% Set MOWOA parameters
params.max_gen = 100;               
params.npop = 100;                 
params.a = 2;                       
params.a_min = 0;                   
params.b = 1;                       

% RUN MOWOA
[final_pop, final_obj, pareto_front] = mowoa(problem, params);

% Draw the Pareto frontier
figure;
plot(pareto_front(:,1), pareto_front(:,2), 'ro', 'MarkerSize', 6);
xlabel('f_1');
ylabel('f_2');
title('MOWOA Pareto frontier');
grid on;
