function [l2] = kernel_l2_no_split(X,Y)
%% Returns the l2 squared distance estimator for the two samples X,Y
%% without data splitting.
%% Uses the bi-variate kernel evaluations as the estimator.

n1 = size(X,2);
n2 = size(Y,2);
d = size(X,1);

%% [est_probs, f, h_old] = kde(X');

%% beta = 4;
%% h = h_old*n1^(-2/(4*beta+d) + 1/(2*beta+d));
beta = d;
h = 0.5*n1^(-2/(4*beta+d) + 1/(2*beta+d));

%% Median heuristic
% $$$ Z = [X Y]';
% $$$ s = size(Z,1);
% $$$ if s > 100
% $$$     Zmed = Z(1:100,:);
% $$$     s = 100;
% $$$ else
% $$$     Zmed = Z;
% $$$ end
% $$$ G = sum((Zmed.*Zmed), 2);
% $$$ Q = repmat(G,1,s);
% $$$ R = repmat(G',s,1);
% $$$ dists = Q+R - 2*Zmed*Zmed';
% $$$ dists = dists - tril(dists);
% $$$ dists = reshape(dists, s^2,1);
% $$$ h = sqrt(0.5*median(dists(dists>0)));

T1 = GaussKernel(h, X');
T2 = GaussKernel(h, Y');
T3 = GaussKernel(h, X', Y');

T1 = 1/(n1*(n1-1)) * sum(sum(T1 - diag(diag(T1))));
T2 = 1/(n2*(n2-1)) * sum(sum(T2 - diag(diag(T2))));
T3 = 2/(n1*n2) * sum(sum(T3));

l2 = T1 + T2 - T3;
