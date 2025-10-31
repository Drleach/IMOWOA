function pop = initialization(npop, nvar, varmin, varmax)
    pop = zeros(npop, nvar);
    for i = 1:npop
        pop(i,:) = varmin + (varmax - varmin) .* rand(1, nvar);
    end
end