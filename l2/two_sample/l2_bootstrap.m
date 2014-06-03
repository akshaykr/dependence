function [h,p,that] = l2_bootstrap(X,Y),
  [that, vals, lb, ub] = bootstrap(X,Y,0.05);
  if 0 >= lb,
      h = 0;
  else,
      h = 1;
  end;

  p = 0;