# IMOWOA
A novel multi-objective whale optimization algorithm (IMOWOA) is proposed. The elite-oriented multi-mode adaptive perturbation (EMAP) and random single-dimensional update (RSDU) strategies are introduced to enable the whale optimization algorithm (WOA) to focus on the elite region and enhance local exploration. Then, a method for calculating the crowding distance is proposed, which integrates multi-scale neighborhoods and nonlinear weights to generate a set of high-quality and evenly distributed non-dominated solutions.
Simply run the main function to start IMOWOA. Pay attention to the difference between crowding_distance1 and crowding_distance2. The calculation formula for the crowding distance in this case is the same as that used in NSGAII. The crowding distance 2, however, is proposed by us.

The new crowding degree formula introduces multi-scale neighborhoods (i±1, i±2), global statistics (μj, σj), and nonlinear terms (exponential terms). While these designs are theoretically more comprehensive, they are also more susceptible to the following factors:

# Parameter Sensitivity: Key parameters are not adapted to your problem.
The new method relies on multiple adjustable parameters, and the default settings may not match your optimization problem, leading to computational bias:
1.Minimum constant ∂., below 1e-5), ∂ may be too large, masking the true differences in τglobal. Conversely, if the target value is large (e.g., above 1e3), ∂ may be too small, causing the denominator (τglobal + ∂) to approach zero and triggering numerical instability (e.g., sudden increases or decreases in crowding degree).
2.Objective weight ω_j: The code defaults to "equal weights," but if the magnitudes of different objectives in your problem vary greatly (e.g., f1 ranges from [0, 1] and f2 ranges from [0, 100]), equal weighting will cause the smaller objectives to be neglected. This results in a crowding degree ranking that favors objectives with larger magnitudes, thus compromising diversity.
3.Boundary handling logic: In the new method, when the neighborhood of i±2 is used with a "small number of individuals on the front" (e.g., only 3 individuals on the front), "neighborhood overlap" occurs (e.g., for i=2, i-2=1 and i+2=3). In such cases, τmid and τlocal become nearly identical, rendering the multi-scale design meaningless and instead increasing computational noise.
# Scenario Adaptability: The assumptions of the new method do not match your problem.
The theoretical advantages of the new method depend on the assumption of "relatively uniform population distribution with no strong correlations between objectives." If your problem does not meet these assumptions, the method’s effectiveness may be compromised:
1.Small population size: The new method relies on global statistics (μj, σj). If the population size is small (e.g., npop&lt;50), μj and σj will be biased due to insufficient sample size, leading to inaccurate computation of τglobal. For example, during the early stages of iteration, when the population’s objective values are concentrated, σj approaches zero, causing the nonlinear term exp(-∆avg/(2σj²)) to approach zero. This results in the crowding degree contribution being "over-suppressed," making it difficult to distinguish individual diversity.
2.Special characteristics of the objective function: If your problem involves strongly correlated multi-objectives (e.g., f2 = 2*f1) or weak objective conflicts, the new method's "multi-scale neighborhood" and "global statistics" may amplify irrelevant differences, leading to the mis-selection of good individuals. For instance, in the ZDT3 problem with a discontinuous Pareto front, the global σ_j in the new method might mask the differences in the local front, making it less precise than the original NSGA-II's "local crowding degree" (which only compares directly neighboring individuals).

# Regarding CEC2020
This is the CEC Multimodal multiobjective optimization (MMO) test problem suite.
----------------------------------------------------------------------------
In the file folder "MMO test problems"
functions: the matlab codes of the test problems
Reference_PSPF_data: the reference PSs and PFs data of the test problems

----------------------------------------------------------------------------

The codes are available on https://github.com/P-N-Suganthan and http://www5.zzu.edu.cn/ecilab/info/1036/1163.htm

Please cite the Technical Report "J. J. Liang, P. N. Suganthan, B. Y. Qu, D. W. Gong, and C. T. Yue, Problem definitions and evaluation criteria for the CEC 2020 special session on multimodal multiobjective optimization, Technical Report, Computational Intelligence Laboratory, Zhengzhou University

----------------------------------------------------------------------------

Any question, please contact jiguang1107@foxmail.com.


