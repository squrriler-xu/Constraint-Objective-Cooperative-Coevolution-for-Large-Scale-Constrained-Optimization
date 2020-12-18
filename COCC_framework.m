function [bestfit, bestphi] = COCC_framework(fname, func_num, Alg, Cht, runs, delta)

% parameters definition
global initial_flag
initial_flag = 0;

global Gen
Gen = 0;

Max_FEs = 3e6;
popsize = 100;
e = 1e-4;

if func_num == 10
    dim = 900;
else
    dim = 1000;
end

[lb, ub] = GetBounds(func_num);

% load groups_fit and groups_phi
[fitGroupsNum, tempFit] = getGroups(func_num, 'groups_fit');
[phiGroupsNum, tempPhi] = getGroups(func_num, 'groups_phi');

path_grouping_result = sprintf('./dg3/result/LSC/f%02d.mat', func_num);
load(path_grouping_result, 'used_FEs');

% initialization population
pop = lb + rand(popsize, dim) * (ub-lb);
[fit, g, h] = feval(fname, pop, func_num);
FEs = popsize + used_FEs;

% 处理约束和约束冲突
g(g < 0) = 0; h = abs(h); h(h <= e) = 0; 
Con = [g, h];
phi = sum(Con, 2);

% initialization epsilon parameter

epsilon = getEpsilon(phi, Max_FEs, Alg);
% CHT = getCHT(Cht, epsilon);

% get best member and best value
P = [pop, fit, phi];
P = sortrows(P, dim+1);
P = sortrows(P, dim+2);
bestfit = P(1, dim+1);
bestphi = P(1, dim+2);
bestmem = P(1, 1:dim);

fprintf('%d, %d, %d\n', func_num, size(tempFit,2), size(tempPhi, 2));

% save the dimension index and contributions for the each groups
% OPTS is parameter of COCC
for i = 1 : fitGroupsNum
    groupsFit(i).index = tempFit{i};       % dimension index
    groupsFit(i).fitDelta = 0;             % contribution
    groupsFit(i).phiDelta = 0;
    groupsFit(i).phi = mean(phi);
    groupsFit(i).OPTS.first = 0;
    groupsFit(i).OPTS.epsilon = epsilon;
end

for i = 1 : phiGroupsNum
    groupsPhi(i).index = tempPhi{i};
    groupsPhi(i).fitDelta = 0;             
    groupsPhi(i).phiDelta = 0;
    groupsPhi(i).phi = mean(phi);
    groupsPhi(i).OPTS.first = 0;
    groupsPhi(i).OPTS.epsilon = epsilon;
end

% parameters of optimizer
iter = 100;

fitFlag = 1;
phiFlag = 1;

fitIndex = 0;
phiIndex = 0;

Flag = '';

CHT.epsilon = epsilon.epsilon0;
CHT.str = Cht;

while (FEs < Max_FEs)
    % 选择器
    [fitFlag, phiFlag, fitIndex, phiIndex, Flag, isChange] =...
        selector(fitFlag, phiFlag, fitIndex, phiIndex, groupsFit, groupsPhi, bestphi, Flag, CHT, delta);
    
    if strcmp(Flag, 'phi')
        dealGroup = groupsPhi(phiIndex);
    elseif strcmp(Flag, 'fit')
        dealGroup = groupsFit(fitIndex);
    end
    
    % 获得子分组
    if isChange == 1
        [gpop, subpop, fit, phi, used_FEs] = getSubPop(fname, func_num, dealGroup, pop, bestmem);
        FEs = FEs + used_FEs;
    end
    
    % 计算剩余代数
    if strcmp(Alg, 'IUDE')
        if (FEs + (iter * popsize * 2) > Max_FEs)
            iter = ceil((Max_FEs - FEs) / (popsize*2));
        end
    else
        if (FEs + (iter * popsize) > Max_FEs)
            iter = ceil((Max_FEs - FEs) / popsize);
        end
    end
    
    % 优化器
    [pop, gpop, subpop, fit, phi, bestmem, bestfit, bestphi, dealGroup, used_FEs, CHT]...
        = optimizer(Alg, func_num, pop, gpop, subpop, fit, phi, bestmem, bestfit, bestphi, lb, ub, iter, dealGroup, CHT);
    FEs = FEs + used_FEs;
    
    if strcmp(Flag, 'phi')
        groupsPhi(phiIndex) = dealGroup;
    elseif strcmp(Flag, 'fit')
        groupsFit(fitIndex) = dealGroup;
    end
    
    fprintf('%s.%s.%2d.%2d| Gen: %d, Func Val: %f, CV: %f\n', Alg, Cht, func_num, runs, Gen, bestfit, bestphi);
end

end

function [fitFlag, phiFlag, fitIndex, phiIndex, Flag, isChange] = selector(fitFlag, phiFlag, fitIndex, phiIndex, groupsFit, groupsPhi, bestphi, Flag, CHT, delta)
% CHT is constraint handing techniques
% isChange 表示是否切换了分组
pf = 0.45;

fitGroupsNum = length(groupsFit);
phiGroupsNum = length(groupsPhi);

fitLastIndex = fitIndex;
phiLastIndex = phiIndex;

lastFlag = Flag;

isChange = 1;

if strcmp(CHT.str, 'EC')
    epsilon = CHT.epsilon;
else
    epsilon = 0;
end

if strcmp(CHT.str, 'EC') || strcmp(CHT.str, 'SF')
    if bestphi <= epsilon
        Flag = 'fit';
    else
        Flag = 'phi';
    end
elseif strcmp(CHT.str, 'SR')
    if bestphi == 0 || rand < pf
        Flag = 'fit';
    else
        Flag = 'phi';
    end
end

if strcmp(Flag, 'fit')
    % 从 fitQ 中选择
    if fitFlag == 1
        fitIndex = fitLastIndex + 1;
        if fitIndex == fitGroupsNum
            fitFlag = 0;
        end
    else
        probFit = (fitGroupsNum : -1 : 1)/sum(1 : fitGroupsNum);
        [~, sortIndex] = sortrows([[groupsFit.fitDelta]', [groupsFit.phiDelta]'], 'descend');
        if bestphi ~= 0
%             i = randi(fitGroupsNum);
%             while rand >= probFit(sortIndex(i))
%                 i = randi(fitGroupsNum);
%             end
            probFit = cumsum(probFit);
            fitIndex = sum(rand > probFit)+1;
        else
            fitIndex = sortIndex(1);
        end
    end
    if fitIndex == fitLastIndex && strcmp(Flag, lastFlag)
        isChange = 0;
    end
else
    % 从 phiQ 中选择
    if phiFlag == 1
        phiIndex = phiLastIndex + 1;
        if phiIndex == phiGroupsNum
            phiFlag = 0;
        end
    else
        if any([groupsPhi.phiDelta] >= delta)
%             while groupsPhi(phiIndex).phiDelta < 1e-4 && rand > 0.1
%                 phiIndex = phiIndex + 1;
%                 if phiIndex > phiGroupsNum
%                     phiIndex = 1;
%                 end
%             end
            geindex = find([groupsPhi.phiDelta] >= delta);
            phiIndex = geindex(randi(length(geindex)));
        else
            phiIndex = randi(phiGroupsNum);
        end
    end
    if phiIndex == phiLastIndex && strcmp(Flag, lastFlag)
        isChange = 0;
    end
end

end

function [pop, gpop, subpop, fit, phi, bestmem, bestfit, bestphi, group, used_FEs, CHT]...
      = optimizer(Alg, func_num, pop, gpop, subpop, fit, phi, bestmem, bestfit, bestphi, lb, ub, max_iter, group, CHT)
  
dim_index  = group.index;

% best value
oldBestfit = bestfit;
oldBestphi = bestphi;

[subpop, fit, phi, newbestsubmem, newbestfit, newbestphi, group.OPTS, used_FEs, epsilon]...
    = feval(Alg, func_num, subpop, gpop, fit, phi, ub, lb, max_iter, dim_index, group.OPTS);

pop(:, dim_index)  = subpop;
gpop(:, dim_index) = subpop;

fitDelta = oldBestfit - newbestfit;
phiDelta = oldBestphi - newbestphi;        % 可能会小于0

if fitDelta < 0
    fitDelta = 0;
end
if phiDelta <= 0
    phiDelta = 0;
end

if (newbestphi < oldBestphi) || (newbestphi == oldBestphi && newbestfit < oldBestfit)
    bestmem(dim_index) = newbestsubmem;
    bestfit = newbestfit;
    bestphi = newbestphi;
end

group.fitDelta = fitDelta;
group.phiDelta = phiDelta;
group.phi = mean(phi);

CHT.epsilon = epsilon;

end

