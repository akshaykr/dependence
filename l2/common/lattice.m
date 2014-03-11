function [lattice] = lattice(d, m)
% Compute all of the points on the integer lattice with
% ||x||_{\infty} \le m
% return a matrix whose rows are the points.
  lattice = zeros((2*m)^d, d);
  s = 1;
  curr = 1;
  while 1
      if curr > (2*m+1)^d;
          return;
      end;
      for i=0:s^d-1
          vec = zeros(1,d);
          for j=1:d
              vec(1,j) = mod(int64(i/s^(j-1)),s);
          end;
          if max(vec) == s-1
              sl = all_sign_pattern(vec);
              lattice(curr:curr+size(sl,1)-1,1:d) = sl;
              curr = curr + size(sl,1);
          end;
      end;
      s = s + 1;
  end;

function [sub_lattice] = all_sign_pattern(v)
% compute all sign patterns for a vector v
  d = size(v,2);
  sub_lattice = zeros(1,d);
  next = 1;
  coords = find(v);
  if isempty(coords),
      return;
  end;
  for i=1:2^(size(coords,2)),
      sgn = zeros(1,size(coords,2));
      for j=1:size(coords,2),
          sgn(1,j) = 2 * mod(int64(i/2^(j-1)), 2) - 1;
      end;
      vec = zeros(1, d);
      for j=1:size(coords,2),
          vec(coords(j)) = sgn(1,j)*v(1,coords(j));
      end;
      sub_lattice(next,1:d) = vec;
      next = next + 1;
  end;