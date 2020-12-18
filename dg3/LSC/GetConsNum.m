function [g_num, h_num, cons_num] = GetConsNum(func_num)
% Get the number of constraint functions 

if func_num == [1, 9]
    g_num = 0;
    h_num = 2;
elseif ismember(func_num, [2, 6, 11])
    g_num = 2;
    h_num = 0;
elseif ismember(func_num, [3, 7])
    g_num = 1;
    h_num = 1;
elseif ismember(func_num, [4])
    g_num = 0;
    h_num = 1;
elseif ismember(func_num, [5])
    g_num = 1;
    h_num = 0;
elseif ismember(func_num, [8, 10, 12])
    g_num = 3;
    h_num = 0;
end

cons_num = g_num + h_num;

end

