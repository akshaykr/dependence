function [l2 lb ub] = confidence_interval(X,Y,alpha, split, ci_method)
%% Compute a confidence interval for the l2 divergence between X and Y. 
%% several options are available:
%% alpha = probability of failure.
%% split = boolean specifying whether to data split or not
%% ci_method = 0 indicates to not data split when estimating the asymptotic variance
               1 indicates to data split
               2 is to use the bootstrap.

  n1 = size(X,2);
  n2 = size(Y,2);
  n = min(n1, n2);

  if split == 1,
      %% data split for estimating L_2
      l2 = kernel_l2(X,Y);
      eff_n = n1/2;
  else,
      %% use all the samples for L_2
      l2 = kernel_l2_no_split(X,Y);
      eff_n = n1;
  end;

  if ci_method == 0,
      %% no data-splitting for variance estimate
      var_est = var_est_no_split(X,Y);
      zalpha = norminv(alpha/2,0,1);

      lb = l2 + sqrt(var_est)*zalpha/sqrt(eff_n);
      ub = l2 - sqrt(var_est)*zalpha/sqrt(eff_n);
  elseif ci_method == 1
      %% data-split variance estimate
      var_est = variance_estimate(X,Y);
      zalpha = norminv(alpha/2,0,1);

      lb = l2 + sqrt(var_est)*zalpha/sqrt(eff_n);
      ub = l2 - sqrt(var_est)*zalpha/sqrt(eff_n);
  elseif ci_method == 2
      %% bootstrap variance estimate
      [theta, vals, lb, ub] = bootstrap(X,Y, alpha);
  end;
