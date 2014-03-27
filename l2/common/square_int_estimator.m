function [that] = square_int_estimator(X),
% Estimates the squared integral of a distribution given iid
% samples from that distribution.

    s = 2;
    n1 = size(X,2);
    d = size(X,1);

    m = 2*int64(max(n1^(2/(4*s+d)), 1));
    Z = double(lattice2(-m:m, d));

    that = 0;
    for i=1:size(Z,1),
        xf = exp(2*pi*1.0j*Z(i,1:d)*X);
        Xf = xf'*xf;
        that = that + 1/(n1*(n1-1))*(sum(sum(Xf)) - sum(diag(Xf)));
    end;