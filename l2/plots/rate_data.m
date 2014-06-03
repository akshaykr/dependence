function [ns, ms, vs] = rate_data(ns, d, iters),
%% Compute the error of the L_2 divergence estimator as a function of
%% the number of samples in d dimension. 
%% returns sample sizes, means and variances of empirical error distribution. 

  theta = 2/(2*sqrt(pi))^d * (1 - exp(-d/4));
  fprintf('theta = %0.4f\n', theta);
  ms = [];
  vs = [];

  for n=ns,
      scores = [];
      for i=1:iters,
          if d == 1,
              x = normrnd(0, 1, [1 n]);
              y = normrnd(1, 1, [1 n]);
          else,
              v = repmat(0, 1, d);
              x = mvnrnd(v, eye(d), n)';
              v = repmat(1, 1, d);
              y = mvnrnd(v, eye(d), n)';
          end;
          [l2] = kernel_l2(x,y);
          scores = [scores abs(l2 - theta)];
      end;
      ms = [ms mean(scores)];
      vs = [vs var(scores)];
      if mod(n, 500) == 0,
          fprintf('n = %d, m = %0.3e, v = %0.3f\n', n, mean(scores), ...
                  var(scores));
      end
  end;
  return;
