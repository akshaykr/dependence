function [h, p, that] = kolmogorov_smirnov(X,Y),
% Kolmogorov Smirnov two-sample test.

d = size(X,1);
if d > 1,
    fprintf('CANNOT RUN KS test in high dimension\n');
    that = -1;
    return;
end;

[h,p,that] = kstest2(X,Y);
return;