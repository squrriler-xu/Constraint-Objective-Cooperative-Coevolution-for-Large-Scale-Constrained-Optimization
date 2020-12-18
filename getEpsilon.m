function epsilon = getEpsilon(phi, Max_FEs, Alg_str)

popsize = length(phi);

if strcmp(Alg_str, 'IUDE')
    TC = ceil(Max_FEs / (popsize*2));
    epsilon0 = max(phi);
    lambda = 6;
    p = 0.5;
    cp = (-log(epsilon0)-lambda) / log(1-p);
    
    epsilon.TC = TC;
    epsilon.epsilon0 = epsilon0;
    epsilon.cp = cp;
    epsilon.p = p;
elseif strcmp(Alg_str, 'DE')
    TC = ceil(Max_FEs / popsize);
    epsilon0 = max(phi);
    lambda = 6;
    p = 0.5;
    cp = (-log(epsilon0)-lambda) / log(1-p);
    
    epsilon.TC = TC;
    epsilon.epsilon0 = epsilon0;
    epsilon.cp = cp;
    epsilon.p = p;
elseif strcmp(Alg_str, 'epsMAgES')
    phi_sort = sort(phi);
    n = ceil(0.9*size(phi_sort, 1)); 
    epsilon0 = mean(phi(1:n));
    TC = ceil(Max_FEs / (popsize*5));
    lambda = 5;
    p = 0.05;
    cp = (-log(epsilon0)-lambda) / log(1-p);
    
    epsilon.TC = TC;
    epsilon.epsilon0 = epsilon0;
    epsilon.cp = cp;
end

end