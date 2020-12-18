clear all; close all;

% 
addpath('LSC');
myfunc = 1 : 12;

func_num = length(myfunc);

%run grouping
% for func=1:length(myfunc)
for func = 1 : 12
    global initial_flag;
    initial_flag = 0;
    %set dimension
    if func == 10
        dim = 900 ;
    else
        dim = 1000;
    end
%     dim = 1000;
    [lb, ub] = GetBounds(func);
    lb = lb * ones(1, dim);
    ub = ub * ones(1, dim);
    
    [groups_fit, FEs] = dg3('LSC', func, lb, ub, dim, 'fit');
    used_FEs = FEs;
    [groups_phi, FEs] = dg3('LSC', func, lb, ub, dim, 'phi');
    used_FEs = used_FEs + FEs;
    filename = sprintf('./result/LSC/f%02d.mat',func);
    
    save(filename,'groups_fit', 'groups_phi', 'used_FEs');

end