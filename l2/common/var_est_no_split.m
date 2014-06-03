function [v] = var_est_no_split(X,Y)
%% Estimate the asymptotic variance of the L_2^2 divergence
%% estimator without splitting the data. This is
%% 4 \Var_{x \sim p}(p(x)) + 4 \Var_{y \sim q}(q(y)) + 
%% 4 \Var_{x \sim p}(q(x)) + 4 \Var_{y \sim q}(p(y))
%% For now we are using a plugin estimator and monte carlo.


  n1 = size(X,2);
  n2 = size(Y,2);
  
  [ep, phat, hp] = kde(X');
  [eq, qhat, hq] = kde(Y');

  T1 = 4*var(phat(X'));
  T2 = 4*var(qhat(Y'));

  T3 = 4*var(qhat(X'));
  T4 = 4*var(phat(Y'));

  v = T1 + T2 + T3 + T4;
