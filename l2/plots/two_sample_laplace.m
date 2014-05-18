addpath('../common/');
addpath('../two_sample/');
addpath('../two_sample/mmd/');
mus = 0.1:0.1:4.0;
[mmdp] = laplace_experiment(mus, @mmd);
[l2p] = laplace_experiment(mus, @l2_test);
[ksp] = laplace_experiment(mus, @kolmogorov_smirnov);


save('two_sample_laplace_exp.mat', 'mus', 'l2p', 'ksp', 'mmdp');
