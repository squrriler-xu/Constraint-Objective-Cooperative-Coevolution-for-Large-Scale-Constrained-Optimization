function [group_num, grouping_result] = getGroups(func_num, str)
% get grouping results from ideal_grouping and combine unidimensional groups

% get grouping results
groups = {};
path_grouping_result = sprintf('./dg3/result/LSC/f%02d.mat', func_num);
load(path_grouping_result, str);

if strcmp('groups_fit', str)
    groups = groups_fit;
else
    groups = groups_phi;
end

% deal grouping results
single_group      = [];
grouping_result   = {};
for i = 1 : size(groups, 2)
    if(length(groups{i}) == 1)
        single_group     = [single_group, groups{i}];
    else
        grouping_result  = [grouping_result, groups(i)];
    end
end

single_num  = length(single_group);
min_size    = 50;
while(single_num > 0)
    cur_size            = min(min_size,single_num);
    grouping_result     = [grouping_result, single_group(1:cur_size)];
    single_group(1 : cur_size) = [];
    single_num          = single_num - cur_size;
end

group_num = size(grouping_result, 2);

end

