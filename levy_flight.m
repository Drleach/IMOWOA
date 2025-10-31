function delta = levy_flight(dim)
    beta = 1.5;
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(beta*gamma((1+beta)/2)*2^((beta-1)/2)))^(1/beta);
    u = randn(1, dim)*sigma;
    v = randn(1, dim);
    delta = u ./ (abs(v).^(1/beta));
end