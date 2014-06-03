
n = 10:10:100;
iters = 200;

fprintf('kernel_l2 no_split bootstrap\n');
[nnb, pnb] = ci_data(ns, 1, iters, 0, 2);

fprintf('kernel_l2 no_split no_split\n');
[nnn, pnn] = ci_data(ns, 1, iters, 0, 0);

fprintf('kernel_l2 split split\n');
[nss, pss] = ci_data(ns, 1, iters, 1, 1);

fprintf('kernel_l2 no_split bootstrap\n');
[nnb, pnb] = ci_data(ns, 1, iters, 0, 2);
