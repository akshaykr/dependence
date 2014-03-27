function [l2] = projDivergence(X,Y)
% Estimates the L_2^2 divergence between distributions generating X
% and Y.

% The estimator is based on three terms all using orthogonal
% series-style estimators.
% for the squared term \int p^2(x) dx we use a U-statistic
% 1/(n*(n-1))\sum_{i \ne j} \sum_{m \in M} \phi_m(X_i) \phi_m(X_j)
% for the cross term we use:
% \frac{1}{n^2}\sum_{i,j} \sum_{m} \phi_m(X_i)\phi_m(Y_j)
% The basis function are the fourier basis

    s = 2;

    n1 = size(X,2);
    n2 = size(Y,2);
    d = size(X,1);

    m = 2*int64(max(min(n1, n2)^(2/(4*s+d)), 1));

    Z = double(lattice2(-m:m, d));

    l2 = 0;
    for i=1:size(Z,1),
        xf = exp(2*pi*1.0j*Z(i,1:d)*X);
        Xf = xf'*xf;
        coeff = 1/(n1*(n1-1))*(sum(sum(Xf)) - sum(diag(Xf)));
% $$$         coeff = 1/(n1*n1)*(sum(sum(Xf)));

        yf = exp(2*pi*1.0j*Z(i,1:d)*Y);
        Yf = yf'*yf;
        coeff = coeff + 1/(n2*(n2-1))*(sum(sum(Yf)) - sum(diag(Yf)));
% $$$         coeff = coeff + 1/(n2*n2)*(sum(sum(Yf)));

        Zf = xf'*yf;
        coeff = coeff - 2/(n1*n2)*(sum(sum(Zf)));
% $$$         fprintf('basis: %d, coeff: %0.4f\n', Z(i,1:d), real(coeff));
        l2 = l2 + real(coeff);
    end;