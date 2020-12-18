function [lbound, ubound] = GetBounds(func_num)
% Get bounds
bounds10 = [12];
bounds50 = [3];
bounds100 = [1, 2, 4, 5, 6, 7, 9, 10, 11];
bounds500 = [8];

if ismember(func_num, bounds10)
    lbound = -10;
    ubound = 10;
elseif ismember(func_num, bounds50)
    lbound = -50;
    ubound = 50;
elseif ismember(func_num, bounds100)
    lbound = -100;
    ubound = 100;
elseif ismember(func_num, bounds500)
    lbound = -500;
    ubound = 500;
end

end

