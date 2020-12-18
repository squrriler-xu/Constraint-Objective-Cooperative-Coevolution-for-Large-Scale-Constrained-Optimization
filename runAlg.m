  
clear all; close all;

run_parallel = 1;

myfunc = 1:12;
func_num = length(myfunc);

if run_parallel == 1
    delete(gcp('nocreate')); 
    parpool('local',func_num);
    spmd(func_num)
        COCC_singleRun(myfunc(labindex));
    end 
else
    %for test
    best = zeros(28, 3);
    for func = 4
        COCC_singleRun(func);
    end
end
