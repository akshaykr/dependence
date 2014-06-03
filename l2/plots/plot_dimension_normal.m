load('dim_two_sample_normal_unscaled.mat');

f = figure();
l1 = plot(ds, l2p, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
l2 = plot(ds, l2perm,  'r-d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
l3 = plot(ds, mmdp, 'g-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
h = xlabel('Dimension (d)');
set(h, 'FontSize', 20);
h = ylabel('Success Probability');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
legend([l1 l2 l3], 'L2asymp', 'L2perm', 'MMD', 4);
saveas(f, './figs/dim_two_sample_unscaled.fig', 'fig');
saveas(f, './figs/dim_two_sample_unscaled.eps', 'epsc');

load('dim_two_sample_normal_scaled.mat');

f = figure();
l1 = plot(ds, l2p, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
l2 = plot(ds, l2perm,  'r-d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
l3 = plot(ds, mmdp, 'g-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
h = xlabel('Dimension (d)');
set(h, 'FontSize', 20);
h = ylabel('Success Probability');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
legend([l1 l2 l3], 'L2asymp', 'L2perm', 'MMD', 4);
saveas(f, './figs/dim_two_sample_scaled.fig', 'fig');
saveas(f, './figs/dim_two_sample_scaled.eps', 'epsc');
