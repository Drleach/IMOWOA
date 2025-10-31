function crowd_dist = crowding_distance(obj, fronts)
% A new crowding degree calculation method based on multi-scale neighborhood and nonlinear weighting
% Input:
%   obj - Matrix of objective values for all individuals (n×num, where n is the number of individuals and num is the number of objectives)
%   fronts - List of non-dominated frontiers (each element is the index of individuals for that frontier)
% Output:
%   crowd_dist - Crowding distance for each individual
    
    [n, num] = size(obj);  % n:个体总数，num:目标函数数量
    crowd_dist = zeros(n, 1);  % 初始化拥挤度
    
    % 1. 计算每个目标函数的全局统计量（均值μ_j和标准差σ_j）
    mu_j = mean(obj);  % μ_j: 第j个目标的均值 (1×num)
    sigma_j = std(obj);  % σ_j: 第j个目标的标准差 (1×num)
    sigma_j(sigma_j == 0) = 1e-6;  % 避免σ_j为0导致后续除零
    
    % 2. 设置目标函数权重（等权重，满足∑ω_j=1）
    omega_j = ones(1, num) / num;  % ω_j: 第j个目标的重要性权重
    
    % 3. 最小常数∂（避免分母为0）
    delta = 1e-6;
    
    % 对每个前沿计算拥挤度
    for f_idx = 1:length(fronts)
        front = fronts{f_idx};  % 当前前沿的个体索引
        front_size = length(front);  % 当前前沿的个体数量
        if front_size == 1
            crowd_dist(front) = Inf;  % 单个个体拥挤度设为无穷大（保持多样性）
            continue;
        end
        
        % 初始化当前前沿的拥挤度
        cd = zeros(front_size, 1);
        
        % 对每个目标函数j计算贡献
        for j = 1:num
            % a. 按第j个目标值对当前前沿个体排序
            [~, sorted_idx] = sort(obj(front, j));  % 升序排序索引
            sorted_front = front(sorted_idx);  % 排序后的个体索引
            f_vals = obj(sorted_front, j);  % 排序后的目标值 (front_size×1)
            
            % b. 计算每个个体在目标j上的多尺度邻域参数
            for k = 1:front_size  % k为排序后的位置（1到front_size）
                % 处理邻域索引（超出范围时取边界值）
                idx_prev1 = max(k-1, 1);  % i-1的位置（最小为1）
                idx_next1 = min(k+1, front_size);  % i+1的位置（最大为front_size）
                idx_prev2 = max(k-2, 1);  % i-2的位置
                idx_next2 = min(k+2, front_size);  % i+2的位置
                
                % 计算τ_local,j(i) = |f_j(i+1) - f_j(i-1)|
                tau_local = abs(f_vals(idx_next1) - f_vals(idx_prev1));
                
                % 计算τ_mid,j(i) = |f_j(i+2) - f_j(i-2)| / 2
                tau_mid = abs(f_vals(idx_next2) - f_vals(idx_prev2)) / 2;
                
                % 计算τ_global,j(i) = |f_j(i) - μ_j| + σ_j
                tau_global = abs(f_vals(k) - mu_j(j)) + sigma_j(j);
                
                % 计算∆_avg,j(i) = (|f_j(i+1)-f_j(i)| + |f_j(i)-f_j(i-1)|) / 2
                delta_avg = (abs(f_vals(idx_next1) - f_vals(k)) + abs(f_vals(k) - f_vals(idx_prev1))) / 2;
                
                % 计算指数项 exp(-∆_avg,j(i)/(2σ_j²))
                exp_term = exp(-delta_avg / (2 * sigma_j(j)^2));
                
                % 计算目标j对个体k的贡献（公式21的项）
                if tau_global + delta < 1e-10  % 避免分母过小
                    term = 0;
                else
                    term = omega_j(j) * (tau_local + tau_mid) / (tau_global + delta) * exp_term;
                end
                
                % 累加至当前个体的拥挤度
                cd(k) = cd(k) + term;
            end
        end
        
        % 将计算好的拥挤度赋值给对应个体（按原前沿索引）
        crowd_dist(sorted_front) = cd;  % 注意sorted_front是排序后的索引，需对应回原个体
    end
end