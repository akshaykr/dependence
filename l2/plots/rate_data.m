function [ns, ms, vs] = rate_data(ns, d, iters),
  theta = 1/(2*sqrt(pi))^d * 2*(exp(1/4) -1)/exp(1/4);
  fprintf('theta = %0.4f\n', theta);
  ms = [];
  vs = [];

  for n=ns,
      scores = [];
      for i=1:iters,
          if d == 1,
              x = normrnd(0, 1, [1 n]);
              y = normrnd(1, 1, [1 n]);
              % x = unifrnd(0, 1, [1, n]);
              % y = unifrnd(0, 1, [1, n]);
              % y(1,:) = y(1,:)/2;
          else,
              v = repmat(0, 1, d);
              x = mvnrnd(v, eye(d), n)';
              v(1) = 1;
              y = mvnrnd(v, eye(d), n)';
              % x = unifrnd(0, 1, [d, n]);
              % y = unifrnd(0, 1, [d, n]);
              % y(1,:) = y(1,:)/2;
          end;
          [l2] = kernel_l2(x,y);
          scores = [scores abs(l2 - theta)];
      end;
      ms = [ms mean(scores)];
      vs = [vs var(scores)];
      fprintf('n = %d, m = %0.3f, v = %0.3f\n', n, mean(scores), var(scores));
  end;
  return;