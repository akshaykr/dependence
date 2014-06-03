%% Script for generating confidence interval data. 
%% Runs the confidence interval experiment for various parameter settings. 

addpath('../common/');

ns = 50:50:2000;
iters = 200;

fprintf('d=1\n');
[ns, ps1] = ci_data(ns, 1, iters, 1, 1);
save('./ci_data_1.mat', 'ns', 'ps1');
fprintf('d=2\n');
[ns, ps2] = ci_data(ns, 2, iters, 1, 1);
save('./ci_data_2.mat', 'ns', 'ps2');
fprintf('d=3\n');
[ns, ps3] = ci_data(ns, 3, iters, 1, 1);
save('./ci_data_3.mat', 'ns', 'ps3');
fprintf('d=5\n');
[ns, ps5] = ci_data(ns, 5, iters, 1, 1);
save('./ci_data_5.mat', 'ns', 'ps5');
fprintf('d=10\n');
[ns, ps10] = ci_data(ns, 10, iters, 1, 1);
save('./ci_data_10.mat', 'ns', 'ps10');


