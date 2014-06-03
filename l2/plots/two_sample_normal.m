addpath('../common/');
addpath('../two_sample/');
addpath('../two_sample/mmd/');
mus = 0.1:0.1:3.0;
d = 3;
l2perm3 = normal_experiment(mus, @l2_permutation, d);
l2p3 = normal_experiment(mus, @l2_test, d);
mmdp3 = normal_experiment(mus, @mmd, d);
l2boot3 = normal_experiment(mus, @l2_bootstrap,d);

%% [ksp] = normal_experiment(mus, @kolmogorov_smirnov);


save('two_sample_normal_exp_d=3.mat', 'mus', 'l2p3', 'l2boot3', 'l2perm3', ...
     'mmdp3');
%%      'ksp', 'mmdp');

d = 10;
l2perm10 = normal_experiment(mus, @l2_permutation, d);
l2p10 = normal_experiment(mus, @l2_test, d);
mmdp10 = normal_experiment(mus, @mmd, d);
l2boot10 = normal_experiment(mus, @l2_bootstrap,d);

%% [ksp] = normal_experiment(mus, @kolmogorov_smirnov);


save('two_sample_normal_exp_d=10.mat', 'mus', 'l2p10', 'l2boot10', 'l2perm10', ...
     'mmdp10');
%%      'ksp', 'mmdp');
