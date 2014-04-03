function [est] = variance_estimate(X,Y),
%% Estimate the asymptotic variance of the L_2^2 divergence
%% estimator. This is
%% 4\Var_{x \sim p}(p(x)) + 4 \Var_{y \sim q}(q(y)) + 
%% \Var_{x \sim p}(q(x)) + \Var_{y \sim q}(p(y))
%% For now we are using a plugin estimator and monte carlo.

  [ep, phat, hp] = kde(X');
  [eq, qhat, hq] = kde(Y');

  T1 = 4*(mean(phat(X').^2) - mean(phat(X'))^2);
  T2 = 4*(mean(qhat(Y').^2) - mean(qhat(Y'))^2);
  T3 = mean(qhat(X').^2) - mean(qhat(X'))^2;
  T4 = mean(phat(Y').^2) - mean(phat(Y'))^2;

  est = T1 + T2 + T3 + T4;