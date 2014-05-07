clear all;
addpath('../two_sample/');
addpath('../common/');
addpath('../two_sample/mmd/');

mus = 0.2:0.2:6.0;
l2p = local_departure(mus, @l2_test);
k2p = local_departure(mus, @kolmogorov_smirnov);
mmdp = local_departure(mus, @mmd);

f = figure();
l1 = plot(mus, l2p, 'blue');
hold on;
l2 = plot(mus, k2p, 'red');
l3 = plot(mus, mmdp, 'green');
legend([l1,l2,l3], 'kde', 'ks', 'mmd', 4);