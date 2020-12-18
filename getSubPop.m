function [gpop, subpop, fit, phi, used_FEs] = GetSubPop(fname, func_num, group, pop, bestmem)
% get subpopulation for each groups
e = 1e-4;

popsize     = size(pop,1);
dim_index   = group.index;
subpop      = pop(:, dim_index);
gpop        = ones(popsize, 1) * bestmem;
gpop(:, dim_index) = subpop;

[fit, g, h] = feval(fname, gpop, func_num);

g(g < 0) = 0; h = abs(h); h(h <= e) = 0; 
Con = [g, h];
phi = sum(Con, 2);

used_FEs    = popsize;
end
