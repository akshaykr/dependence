% Unit Test for genData.m

% Summary
fprintf('p1: %f, sum(mixture_idxs) = %d\n', p1, sum(mixture_idxs));
fprintf('MU1: %s, mean(mix1): %s\n', mat2str(MU1), mat2str(mean(mix1_means)));
fprintf('MU2: %s, mean(mix2): %s\n', mat2str(MU2), mat2str(mean(mix2_means)));

% Plot each distribution
for i = 1:NUM_DISTS
  curr_data = Data{i};
  if NUM_DIMS == 1, hist(curr_data);
  else, plot(curr_data(:,1), curr_data(:,2), 'x');
  end
  axis([-2 12, -2, 12]);
  title_str = sprintf('p1: %f, dims: %d, mix_idx: %d, mean: (%.2f, %.2f),', ...
    p1, NUM_DIMS, mixture_idxs(i), dist_means(i,1), dist_means(i,2) );
  title(title_str);
  pause;
end

