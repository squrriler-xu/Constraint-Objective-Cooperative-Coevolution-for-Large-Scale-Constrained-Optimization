 function COCC_singleRun(func)

runs = 25;
% run
% bestval = zeros(1, runs);
% bestPhi = zeros(1, runs);

% strategy = 'DE';
% strategy = 'IUDE';
% strategy = 'epsMAgES';
fname = 'LSC';
% Cht = 'SF';
% Cht = 'EC';
% Cht = 'SR';

for a = 1:3
    if a == 1
        strategy = 'DE';
    elseif a == 2
        strategy = 'IUDE';
    elseif a == 3
        strategy = 'epsMAgES';
    end
    
    for b = 1:3
        if b == 1
            Cht = 'SF';
        elseif b == 2
            Cht = 'EC';
        elseif b == 3
            Cht = 'SR';
        end
%         for delta = [1e-1, 1e-2, 1e-3, 1e-4]

            delta = 1e-1;
            bestval = zeros(1, runs);
            bestPhi = zeros(1, runs);
            for i = 1:runs
                s1 = RandStream('mt19937ar', 'Seed', i);
                RandStream.setGlobalStream(s1);
                [bestval(i), bestPhi(i)] = COCC_framework(fname, func, strategy, Cht, i, delta);
                %     [bestval(i), bestPhi(i)] = eCCPSO_DG(fname, func);
                fprintf('%2d.%2d| Func Val: %f, CV: %f\n', func, i, bestval(i), bestPhi(i));
            end
            
            mean_val = mean(bestval);
            mean_Phi = mean(bestPhi);
            
            fr = sum(bestPhi == 0) / runs;
            
            filename=sprintf('./result/COCC_%s_%s_f%02d.csv', strategy, Cht, func);
            
            % filename=sprintf('./result/eCCPSO_DG_f%02d.csv', func);
            output=[func, mean_val, mean_Phi, fr];
            csvwrite(filename,output);
    end
end
 end