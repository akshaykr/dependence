function [h,p,that] = l2_bootstrap(X,Y),
%% bootstrapped two-sample test for L2 divergence.

  [that, vals, lb, ub] = bootstrap(X,Y,0.05);
  if 0 >= lb,
      h = 0;
  else,
      h = 1;
  end;

  p = 0;
