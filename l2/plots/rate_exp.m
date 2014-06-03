%% Script for rate of convergence experiments. 

addpath('../common/');

ns = 50:50:10000;
iters = 10;

fprintf('d = 1\n');
[ns, ms1, vs1] = rate_data(ns, 1, iters);
save('./rate_data_1.mat', 'ns', 'ms1', 'vs1');
fprintf('d = 2\n');
[ns, ms2, vs2] = rate_data(ns, 2, iters);
save('./rate_data_2.mat', 'ns', 'ms2', 'vs2');
fprintf('d = 3\n');
[ns, ms3, vs3] = rate_data(ns, 3, iters);
save('./rate_data_3.mat', 'ns', 'ms3', 'vs3');
fprintf('d = 5\n');
[ns, ms5, vs5] = rate_data(ns, 5, iters);
save('./rate_data_5.mat', 'ns', 'ms5', 'vs5');
fprintf('d = 10\n');
[ns, ms10, vs10] = rate_data(ns, 10, iters);
save('./rate_data_10.mat', 'ns', 'ms10', 'vs10');

save('./rate_data.mat', 'ns', 'ms1', 'vs1', 'ms2', 'vs2', 'ms3', ...
     'vs3', 'ms5', 'vs5', 'ms10', 'vs10');
