function [v] = var_est_no_split(X,Y)

  n1 = size(X,2);
  n2 = size(Y,2);
  
  [ep, phat, hp] = kde(X');
  [eq, qhat, hq] = kde(Y');

  T1 = 4*var(phat(X'));
  T2 = 4*var(qhat(Y'));

  T3 = 4*var(qhat(X'));
  T4 = 4*var(phat(Y'));

  v = T1 + T2 + T3 + T4;