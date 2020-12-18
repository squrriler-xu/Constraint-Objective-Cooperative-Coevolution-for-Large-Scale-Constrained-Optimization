% The SaNSDE algorithm can be found in:
% Zhenyu Yang, Ke Tang and Xin Yao, "Self-adaptive Differential Evolution with
% Neightborhood Search", in Proceedings of the 2008 IEEE Congress on 
% Evolutionary Computation (CEC2008), Hongkong, China, 2008, pp. 1110-1116.

function [pop, fit, phi, bestmem, bestfit, bestphi, OPTS, FEs, epsilon]...
    = DE(func_num, pop, gpop, fit, phi, ub, lb, max_iter, dim_index, OPTS)

fname = 'LSC';
e = 1e-4;

global Gen
iter = 0;

[ps, D] = size(pop);

FEs = 0;

epsilon0 = OPTS.epsilon.epsilon0;
cp = OPTS.epsilon.cp;
Max_gen = OPTS.epsilon.TC;
pc = OPTS.epsilon.p;

epsilon = 0;

while iter < max_iter
iter = iter + 1;
Gen = Gen + 1;

if Gen/Max_gen <= pc
    epsilon = epsilon0*((1 - Gen/Max_gen)^cp);
else
    epsilon = 0;
end

% rot = (0:1:ps-1);               % rotating index array (size NP)
% 
% % 索引指针数组
% ind = randperm(4);              % index pointer array
% 
% % 打乱向量位置，并根据ind进行旋转
% a1  = randperm(ps);             % shuffle locations of vectors
% rt = rem(rot+ind(1),ps);        % rotate indices by ind(1) positions
% a2  = a1(rt+1);                 % rotate vector locations
% rt = rem(rot+ind(2),ps);
% a3  = a2(rt+1);

index = zeros(ps, 3);
for i = 1:ps
    index(i, :) = randperm(ps-1, 3);
    index(i, index(i, :) >= i) = index(i, index(i, :) >= i) + 1;
end

% 打乱种群
pm1 = pop(index(:, 1),:);             % shuffled population 1
pm2 = pop(index(:, 2),:);             % shuffled population 2
pm3 = pop(index(:, 3),:);             % shuffled population 3

%% 更新 CR
CR = 0.9;

%% 更新 F
F = 0.5;

%% 变异交叉
aa = rand(ps, D) < CR;
index = find(sum(aa') == 0);
tmpsize = size(index, 2);
for k=1:tmpsize
    bb = ceil(D*rand);
    aa(index(k), bb) = 1;
end
        
mui = aa;
mpo = mui < 0.5;                % inverse mask to mui

ui = pm3 + F .* (pm1 - pm2);
ui = pop .* mpo + ui .* mui;


%% 越界处理
reflect = find(ui > ub);
ui(reflect) = ub - mod((ui(reflect) - ub), (ub - lb));

reflect = find(ui < lb);
ui(reflect) = lb + mod((lb - ui(reflect)), (ub - lb));

gpop(:, dim_index) = ui;

[tempfit, g, h] = feval(fname, gpop, func_num);     % 评估种群适应值
g(g < 0) = 0; h(abs(h) <= e) = 0;
tempCon = [g, h]; h = abs(h);
tempphi = sum([g, h], 2);
FEs = FEs + ps;

Phi_improve = (tempphi > epsilon | phi > epsilon) & tempphi < phi;
fit_improve = (tempphi <= epsilon & phi <= epsilon | (tempphi == phi)) & tempfit <= fit;
improve = Phi_improve | fit_improve;

pop(improve, :) = ui(improve, :);
fit(improve) = tempfit(improve);
phi(improve) = tempphi(improve);

end
P = [pop, fit, phi;];
P = sortrows(P, D+1);
P = sortrows(P, D+2);
bestfit = P(1, D+1);
bestphi = P(1, D+2);
bestmem = P(1, 1:D);

end
