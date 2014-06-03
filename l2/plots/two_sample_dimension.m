addpath('../common/');
addpath('../two_sample/');
addpath('../two_sample/mmd/');

mu = 1.0;
ds = 5:5:100;

l2p = dim_normal_experiment(mu, @l2_test, ds, 1);
l2perm = dim_normal_experiment(mu, @l2_permutation, ds, 1);
mmdp = dim_normal_experiment(mu, @mmd, ds, 1);

save('./dim_two_sample_normal_scaled.mat', 'ds', 'l2perm', 'l2p', 'mmdp');

l2p = dim_normal_experiment(mu, @l2_test, ds, 0);
l2perm = dim_normal_experiment(mu, @l2_permutation, ds, 0);
mmdp = dim_normal_experiment(mu, @mmd, ds, 0);

save('./dim_two_sample_normal_unscaled.mat', 'ds', 'l2perm', 'l2p', 'mmdp');
