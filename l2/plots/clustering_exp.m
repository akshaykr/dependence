addpath('../common/');

n = 25;
ms = 10:10:100;
mus = 0.25:0.25:1.0;
d = 5;
errs = zeros(size(mus, 2), size(ms,2));
iters = 5;
for j=1:size(mus,2),
    for k = 1:size(ms,2),
        err = [];
        for i=1:iters,
            [Data, Labels] = mixture_model(n,ms(k),d,mus(j),0.1);
            [c1, c2] = cluster_mixture(Data);
            s1 = size(setdiff(c1, find(Labels == 1)),1);
            s2 = size(setdiff(c1, find(Labels == -1)),1);
            err = [err min(s1, s2)/n];
        end
        fprintf('mu = %0.2f m = %d err = %0.2f\n', mus(j), ms(k), mean(err));
        errs(j,k) = mean(err);
    end;
end;