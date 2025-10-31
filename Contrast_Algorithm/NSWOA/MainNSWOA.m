% NSWOA
clc
clear all
D = 30;
M = 2;
K=M+D;
LB = ones(1, D).*0;
UB = ones(1, D).*1; 
Max_iteration = 100;  
SearchAgents_no = 100;     
ishow = 10;
%% Initialize the population
chromosome = initialize_variables(SearchAgents_no, M, D, LB, UB);

%% Sort the initialized population
intermediate_chromosome = non_domination_sort_mod(chromosome, M, D);

%% Perform Selection
Population = replace_chromosome(intermediate_chromosome, M,D,SearchAgents_no);

%% Start the evolution process
Pareto = NSWOA(D,M,LB,UB,Population,SearchAgents_no,Max_iteration,ishow);
save Pareto.txt Pareto -ascii;  % save data for future use

%% Plot data
if M == 2
    plot_data2(M,D,Pareto)
elseif M == 3
    plot_data_TCQ(M,D,Pareto); 
end