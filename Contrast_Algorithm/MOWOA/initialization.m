% The following auxiliary functions are the same as those implemented in NSGA-II and can be directly reused.
function pop = initialization(npop, nvar, varmin, varmax)
    pop = zeros(npop, nvar);
    for i = 1:npop
        pop(i,:) = varmin + (varmax - varmin) .* rand(1, nvar);
    end
end