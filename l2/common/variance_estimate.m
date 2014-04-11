function [est] = variance_estimate(X,Y),
%% Estimate the asymptotic variance of the L_2^2 divergence
%% estimator. This is
%% 4\Var_{x \sim p}(p(x)) + 4 \Var_{y \sim q}(q(y)) + 
%% \Var_{x \sim p}(q(x)) + \Var_{y \sim q}(p(y))
%% For now we are using a plugin estimator and monte carlo.

  [ep, phat, hp] = kde(X');
  [eq, qhat, hq] = kde(Y');
  
  n1 = size(X,2);
  n2 = size(Y,2);

  X1 = X(:,1:n1/2);
  X2 = X(:,n1/2:end);
  Y1 = Y(:,1:n2/2);
  Y2 = Y(:,n2/2:end);

  %% T1 = 4*(mean(phat(X1').^2) - mean(phat(X1'))^2);
  %% T2 = 4*(mean(qhat(Y1').^2) - mean(qhat(Y1'))^2);
  %% T3 = mean(qhat(X2').^2) - mean(qhat(X2'))^2;
  %% T4 = mean(phat(Y2').^2) - mean(phat(Y2'))^2;
  T1 = 4*var(phat(X1'));
  T2 = 4*var(qhat(Y1'));
  T3 = 4*var(qhat(X2'));
  T4 = 4*var(phat(Y2'));

  est = T1 + T2 + T3 + T4;