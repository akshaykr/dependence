function [l2 lb ub] = confidence_interval(X,Y,alpha)

  n1 = size(X,2);
  n2 = size(Y,2);

  n = min(n1, n2);

  l2 = kernel_l2(X,Y);
  var_est = variance_estimate(X,Y);

  zalpha = norminv(alpha/2,0,1);

  lb = l2 + sqrt(var_est)*zalpha/sqrt(n1/2);
  ub = l2 - sqrt(var_est)*zalpha/sqrt(n1/2);
