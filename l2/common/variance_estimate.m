function [est] = variance_estimate(X,Y),
%% Estimate the asymptotic variance of the L_2^2 divergence
%% estimator. This is
%% 4 \Var_{x \sim p}(p(x)) + 4 \Var_{y \sim q}(q(y)) + 
%% 4 \Var_{x \sim p}(q(x)) + 4 \Var_{y \sim q}(p(y))
%% For now we are using a plugin estimator and monte carlo.

  n1 = size(X,2);
  n2 = size(Y,2);
  X1 = X(:,1:n1/2);
  X2 = X(:,n1/2:end);
  Y1 = Y(:,1:n2/2);
  Y2 = Y(:,n2/2:end);

  [ep, phat1, hp] = kde(X1');
  [eq, qhat1, hq] = kde(Y1');
  T1 = 4*var(phat1(X1'));
  T2 = 4*var(qhat1(Y1'));

  [ep, phat2, hp] = kde(X2');
  [eq, qhat2, hq] = kde(Y2');
  T3 = 4*var(qhat2(X2'));
  T4 = 4*var(phat2(Y2'));

  est = T1 + T2 + T3 + T4;