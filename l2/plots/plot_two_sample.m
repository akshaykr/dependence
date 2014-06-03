%% Script for plotting two sample data (as a function of mean separation)

load('two_sample_normal_exp.mat');

f = figure();
l1 = plot(mus, l2p, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
l2 = plot(mus, ksp,  'r-d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
l3 = plot(mus, mmdp, 'g-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
l4 = plot(mus, l2boot, 'k-d', 'LineWidth', 2, 'MarkerSize', 8, ...
          'MarkerFaceColor', 'k');
l5 = plot(mus, l2perm, 'm-o', 'LineWidth', 2, 'MarkerSize', 8, ...
          'MarkerFaceColor', 'm');

h = xlabel('Mean Separation (mu)');
set(h, 'FontSize', 20);
h = ylabel('Success Probability');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
legend([l1 l2 l3 l4 l5], 'L2asymp', 'KS', 'MMD', 'L2boot', 'L2perm', 4);
saveas(f, './figs/two_sample_normal.fig', 'fig');
saveas(f, './figs/two_sample_normal.eps', 'epsc');

load('two_sample_normal_exp_d=3.mat');

f = figure();
l1 = plot(mus, l2p3, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
l3 = plot(mus, mmdp3, 'g-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
l4 = plot(mus, l2boot3, 'k-d', 'LineWidth', 2, 'MarkerSize', 8, ...
          'MarkerFaceColor', 'k');
l5 = plot(mus, l2perm3, 'm-o', 'LineWidth', 2, 'MarkerSize', 8, ...
          'MarkerFaceColor', 'm');

h = xlabel('Mean Separation (mu)');
set(h, 'FontSize', 20);
h = ylabel('Success Probability');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
legend([l1 l3 l4 l5], 'L2asymp', 'MMD', 'L2boot', 'L2perm', 4);
saveas(f, './figs/two_sample_normal_d=3.fig', 'fig');
saveas(f, './figs/two_sample_normal_d=3.eps', 'epsc');

load('two_sample_normal_exp_d=10.mat');

f = figure();
l1 = plot(mus, l2p10, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
l3 = plot(mus, mmdp10, 'g-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
l4 = plot(mus, l2boot10, 'k-d', 'LineWidth', 2, 'MarkerSize', 8, ...
          'MarkerFaceColor', 'k');
l5 = plot(mus, l2perm10, 'm-o', 'LineWidth', 2, 'MarkerSize', 8, ...
          'MarkerFaceColor', 'm');

h = xlabel('Mean Separation (mu)');
set(h, 'FontSize', 20);
h = ylabel('Success Probability');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
legend([l1 l3 l4 l5], 'L2asymp', 'MMD', 'L2boot', 'L2perm', 2);
saveas(f, './figs/two_sample_normal_d=10.fig', 'fig');
saveas(f, './figs/two_sample_normal_d=10.eps', 'epsc');
