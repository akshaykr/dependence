addpath('../common/');
addpath('../two_sample/');
addpath('../two_sample/mmd/');
mus = 0.1:0.1:4.0;
[mmdp] = normal_experiment(mus, @mmd);
[l2p] = normal_experiment(mus, @l2_test);
[ksp] = normal_experiment(mus, @kolmogorov_smirnov);


save('two_sample_normal_exp.mat', 'mus', 'l2p', 'ksp', 'mmdp');
