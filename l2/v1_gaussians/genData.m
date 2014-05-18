% Constants
DIRICHLET_PRIOR = [2; 2];

% Define Parameters for the prior distributions
MU1 = 0 * ones(NUM_DIMS, 1);
S1 = eye(NUM_DIMS);
MU2 = 10 * ones(NUM_DIMS, 1);
S2 = eye(NUM_DIMS);
% Covariance matrix of each distribution
SIGMA = eye(NUM_DIMS);

Data = {};

% Sample from mixture
% mixture_probs = dirichlet_sample(DIRICHLET_PRIOR, 1);
% p1 = mixture_probs(1);
% p2 = mixture_probs(2);
p1 = 0.6;
p2 = 0.4;
mixture_idxs = double(rand(NUM_DISTS, 1) < p1);
% Generate NUM_DISTS random means from each prior
mix1_means = bsxfun(@plus, randn(NUM_DISTS, NUM_DIMS) * chol(S1), MU1');
mix2_means = bsxfun(@plus, randn(NUM_DISTS, NUM_DIMS) * chol(S2), MU2');
% Choose which mean to pick depending on mixture_idxs
dist_means = bsxfun(@times, mix1_means, mixture_idxs) + ...
             bsxfun(@times, mix2_means, (1 - mixture_idxs) );

% Now sample from each distribution
for dist_iter = 1:NUM_DISTS
  Data{dist_iter} = bsxfun(@plus, ...
                           randn(NUM_PTS, NUM_DIMS) * chol(SIGMA), ...
                           dist_means(dist_iter) );
end

