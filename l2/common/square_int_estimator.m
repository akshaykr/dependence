function [that] = square_int_estimator(X),
% Estimates the squared integral of a distribution given iid
% samples from that distribution.

    n1 = size(X,2);
    d = size(X,1);

    beta = d;

    h = n1^(-2/(4*beta+d));
    
    T1 = GaussKernel(h, X');

    that = 1/(n1*(n1-1))*sum(sum(T1 - diag(diag(T1))));