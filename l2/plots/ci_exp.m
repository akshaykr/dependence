addpath('../common/');

ns = 12:4:100;
iters = 50;

fprintf('d=1\n');
[ns, ps1] = ci_data(ns, 1, iters);
fprintf('d=2\n');
[ns, ps2] = ci_data(ns, 2, iters);
fprintf('d=3\n');
[ns, ps3] = ci_data(ns, 3, iters);
save('./c_data_1_2_3.mat', 'ns', 'ps1', 'ps2', 'ps3');

