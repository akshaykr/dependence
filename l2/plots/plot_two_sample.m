load('two_sample_laplace_exp.mat');

f = figure();
l1 = plot(mus, l2p, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
l2 = plot(mus, ksp,  'r-d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
l3 = plot(mus, mmdp, 'g-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');

h = xlabel('Mean Separation (mu)');
set(h, 'FontSize', 20);
h = ylabel('Success Probability');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
legend([l1 l2 l3], 'KDE', 'KS', 'MMD', 4);
saveas(f, './figs/two_sample_laplace.fig', 'fig');
saveas(f, './figs/two_sample_laplace.eps', 'epsc');
