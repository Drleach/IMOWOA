function f = NSWOA(D,M,LB,UB,Pop,SearchAgents_no,Max_iteration,ishow)
%% IMOWOA
%% Algorithm Variables
K = D+M;
Whale_pos = Pop(:,1:K+1);
Whale_pos_ad = zeros(SearchAgents_no,K);
%% Optimization Circle
Iteration = 1;
while Iteration<=Max_iteration % for each generation    
    for i = 1:SearchAgents_no   %  (Moth for each individual)
        j = floor(rand()* SearchAgents_no) + 1;
        while j==i
            j = floor(rand()* SearchAgents_no) + 1;
        end
        SF=round(1+rand); %% Scaling factor to perform best coverage in MFO
        % randomly select the best organism from first non-dominated front of Whale_pos
        Ffront = Whale_pos((find(Whale_pos(:,K+1)==1)),:); % first front
        ri =  floor(size(Ffront,1)*rand())+1; % ri is random index
        sorted_population = Whale_pos(ri,1:D);
        % Calculate new solution 
        Whale_posNew1 = Whale_pos(i,1:D)+rand(1,D).*(sorted_population-SF.*Whale_pos(i,1:D));
        % Handling constraints
        Whale_posNew1 = bound(Whale_posNew1(:,1:D),UB,LB); 
        % Evaluate function at Whale_posNew1 
        Whale_posNew1(:,D + 1: K) = evaluate_objective(Whale_posNew1(:,1:D));
        % For the first trail Whale_posNew1
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (Whale_posNew1(:,D+k)<Whale_pos(i,D+k))
                dom_less = dom_less + 1;
            elseif (Whale_posNew1(:,D+k)== Whale_pos(i,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Whale_posNew1 dominates
            % target vector Xi. Replace Xi by Whale_posNew1 in current Whale_possystem and
            % add Xi to advanced population 2
            Whale_pos_ad(i,1:K) = Whale_pos(i,1:K); % add Xi to advanced Whale_pos
            Whale_pos(i,1:K) = Whale_posNew1(:,1:K); % replace Xi by Whale_posNew1            
        else % else Add Xi (trial vector) to advanced Whale_pos
            Whale_pos_ad(i,1:K)= Whale_posNew1;
        end % end if
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
                dom_more = dom_more +1;           
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Whale_posNew1 dominates
            Whale_pos_ad(j,1:K) = Whale_pos(j,1:K); % add Xi to advanced Whale_pos
        end % end if
        j = floor(rand()* SearchAgents_no) + 1;
        while j==i
            j = floor(rand()* SearchAgents_no) + 1;
        end 
        a=2-Iteration*((2)/Max_iteration ); % a decreases linearly fron 2 to 0 in Eq. (2.3)
        a2=-1+Iteration*((-1)/Max_iteration );
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]       
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper 
        b=1;               %  parameters in Eq. (2.5)
        t=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        p = rand();        % p in Eq. (2.6) 
            if p<0.5 % Update the position of the moth with respect to its corresponsing flame
        % Calculate new solution       
        X_rand = sorted_population;
        Whale_posNew1 = Whale_pos(i,1:D)+X_rand-A.*abs(C*X_rand-Whale_pos(j,1:D));
            elseif p>=0.5  
        Whale_posNew1 = Whale_pos(i,1:D)+abs(sorted_population-Whale_pos(j,1:D))*exp(b.*t).*cos(t.*2*pi)+sorted_population;
            end
        Whale_posNew1 = bound(Whale_posNew1(:,1:D),UB,LB);
        Whale_posNew1(:,D + 1: K) = evaluate_objective(Whale_posNew1(:,1:D));
        % Nondomination checking of trial individual
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (Whale_posNew1(:,D+k)<Whale_pos(i,D+k))
                dom_less = dom_less + 1;
            elseif (Whale_posNew1(:,D+k)== Whale_pos(i,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Whale_posNew1 dominates
            % target vector Xi. Replace Xi by Whale_posNew1 in current Whale_possystem and
            % add Xi to advanced population
            Whale_pos_ad(i,1:K) = Whale_pos(i,1:K); % add Xi to advanced Whale_pos
            Whale_pos(i,1:K) = Whale_posNew1(:,1:K); % replace Xi by Whale_posNew1            
        else % else Add Xi (trial vector) to advanced Whale_pos
            Whale_pos_ad(i,1:K)= Whale_posNew1;
        end % end if
        j = floor(rand()* SearchAgents_no) + 1;
        while j==i
            j = floor(rand()* SearchAgents_no) + 1;
        end
        parasiteVector=Whale_pos(i,1:D);
        seed=randperm(D);
        pick=seed(1:ceil(rand*D));  % select random dimension
        parasiteVector(:,pick)=rand(1,length(pick)).*(UB(pick)-LB(pick))+LB(pick);        
        % Evaluate the Parasite Vector
        parasiteVector(:,D + 1: K) = evaluate_objective(parasiteVector(:,1:D));
        % Nondomination checking of trial individual
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1:M
            if (parasiteVector(:,D+k)<Whale_pos(j,D+k))
                dom_less = dom_less + 1;
            elseif (parasiteVector(:,D+k)== Whale_pos(j,D+k))
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more +1;
            end
        end % end for k
        if dom_more == 0 && dom_equal ~= M %  If trial vector Whale_posNew1 dominates
            % target vector Xi. Replace Xi by Whale_posNew1 in current Whale_possystem and
            % add Xi to advanced population
            Whale_pos_ad(j,1:K) = Whale_pos(j,1:K); % add Xi to advanced Whale_pos
            Whale_pos(j,1:K) = parasiteVector(:,1:K); % replace Xi by Whale_posNew1            
        else % else Add Xi (trial vector) to advanced Whale_pos
            Whale_pos_ad(j,1:K)= parasiteVector;
        end % end if 
    end % end for i
    if rem(Iteration, ishow) == 0
        fprintf('Generation: %d\n', Iteration);        
    end
    Whale_pos_com = [Whale_pos(:,1:K) ; Whale_pos_ad];
    intermediate_Whale_pos = non_domination_sort_mod(Whale_pos_com, M, D);
    Pop  = replace_chromosome(intermediate_Whale_pos, M,D,SearchAgents_no);
    Whale_pos=Pop(:,1:K+1); %    
    Iteration = Iteration+1;
end 
f= Whale_pos;

% Check the boundary limit
function a=bound(a,ub,lb)
a(a>ub)=ub(a>ub); a(a<lb)=lb(a<lb);