function ShiftedAndRotated
% checked
% p 打乱维度
% s 分组集合数量
% M_10, M_30, M_50, M_100 不同分组对应的旋转矩阵
% w 不同分组对应的权重
% o 位移

Ast = RandStream('mt19937ar', 'Seed', 3);
RandStream.setGlobalStream(Ast);

for func_num = 1 : 12
    if func_num == 10
        D = 900;
    else
        D = 1000;
    end
    [lb, ub] = GetBounds(func_num);
    o = computeShiftVector(D, lb, ub);
    
    p = randperm(D);
    
    M_10 = computeRotation(10);
    M_30 = computeRotation(30);
    M_50 = computeRotation(50);
    M_100 = computeRotation(100);
    M = computeRotation(1000);
    
    filename = sprintf('Function%02d.mat', func_num);
    % 1. Fully-separable Functions
    if     (func_num ==  1)
        save(filename, 'o');
    % 2. Partially Additively Separable Functions
    %    2.1. Functions with a separable subcomponent:
	elseif (func_num ==  2)
        s = [30, 10, 30, 50, 50, 30, 100];
        w = 10.^(2*randn(1, length(s)));
        save(filename, 'o', 'M_10', 'M_30', 'M_50', 'M_100', 's', 'p', 'w');
    elseif (func_num ==  3)
        s = [50, 100, 10, 30, 10, 100, 30, 50, 30, 100, 10, 100, 50, 50, 30, 10, 30, 50, 100, 10, 50];
        w = 10.^(2*randn(1, length(s)));
        save(filename, 'o', 'M_10', 'M_30', 'M_50', 'M_100', 's', 'p', 'w');
    %    2.2. Functions with no separable subcomponents:
    elseif (func_num ==  4)
        s = [30, 50, 100, 10, 30, 10, 100, 50, 30, 100, 10, 100, 50, 50, 30, 10, 30, 50, 100, 10, 50];
        w = 10.^(2*randn(1, length(s)));
        save(filename, 'o', 's', 'p', 'w');
	elseif (func_num ==  5)
        s = [30, 50, 100, 10, 30, 10, 100, 50, 30, 100, 10, 100, 50, 50, 30, 10, 30, 50, 100, 10, 50];
        w = 10.^(2*randn(1, length(s)));
        save(filename, 'o', 'M_10', 'M_30', 'M_50', 'M_100', 's', 'p', 'w');
    elseif (func_num ==  6)
        M2_10 = computeRotation(10);
        M2_30 = computeRotation(30);
        M2_50 = computeRotation(50);
        M2_100 = computeRotation(100);
        s = [30, 50, 100, 10, 30, 10, 100, 50, 30, 100, 10, 100, 50, 50, 30, 10, 30, 50, 100, 10, 50];
        w = 10.^(2*randn(1, length(s)));
        save(filename, 'o', 'M_10', 'M_30', 'M_50', 'M_100', 'M2_10', 'M2_30', 'M2_50', 'M2_100', 's', 'p', 'w');
    elseif (func_num ==  7)
        s = [30, 50, 100, 10, 30, 10, 100, 50, 30, 100, 10, 100, 50, 50, 30, 10, 30, 50, 100, 10, 50];
        w = 10.^(2*randn(1, length(s)));
        save(filename, 'o', 'M_10', 'M_30', 'M_50', 'M_100', 's', 'p', 'w');
    elseif (func_num ==  8)
        s = [30, 50, 100, 10, 30, 10, 100, 50, 30, 100, 10, 100, 50, 50, 30, 10, 30, 50, 100, 10, 50];
        w = 10.^(2*randn(1, length(s)));
        save(filename, 'o', 'M_10', 'M_30', 'M_50', 'M_100', 's', 'p', 'w');
%         %3. Overlapping Functions
    elseif (func_num == 9)
        save(filename, 'o');
    elseif (func_num == 10)
        m = 5;
        s = [30, 50, 100, 10, 30, 10, 100, 50, 30, 100, 10, 100, 50, 50, 30, 10, 30, 50, 100, 10, 50];
        save(filename, 'o', 'm', 'M_10', 'M_30', 'M_50', 'M_100', 's', 'p');
%         % 4. Fully Non-separable Functions
    elseif (func_num == 11)
        s = [30, 50, 100, 10, 30, 10, 100, 50, 30, 100, 10, 100, 50, 50, 30, 10, 30, 50, 100, 10, 50];
        save(filename, 'o', 'M_10', 'M_30', 'M_50', 'M_100', 'M', 's', 'p');
    elseif (func_num == 12)
        save(filename, 'o', 'M');
    end
end
end

function [Q R] = computeRotation(D)
    A = randn(D);

    Q = zeros(D);
    R = zeros(D);
    for j = 1:D
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end

end

function s = computeShiftVector(D, lb, ub)
    s = zeros(D, 1);
    hw = (ub - lb) / 2.0;
    middle = lb + hw;
    for i=1:D
        s(i) = middle + randn * hw;
        while((s(i) < lb) || (s(i) > ub))
            s(i) = middle + randn * hw;
        end
    end
    s = s';
end

