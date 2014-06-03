%% load('./rate_data.mat');
load('./rate_data_1.mat');
load('./rate_data_3.mat');
load('./rate_data_5.mat');
load('./rate_data_10.mat');

r = 2*sqrt(pi); %% l_2^2 of 1-d gaussian density.

f = figure();
l1 = errorbar(ns, ms1*r, vs1*r^2, 'b-o', 'MarkerSize', 4, ...
              'MarkerFaceColor', 'b', 'LineWidth', 1);
hold on;
l2 = errorbar(ns, ms3*r^3, vs3*r^6, 'k-d', 'MarkerSize', 4, ...
              'MarkerFaceColor', 'k', 'LineWidth', 1);
l3 = errorbar(ns, ms5*r^5, vs5*r^10, 'r-d', 'MarkerSize', 4, ...
              'MarkerFaceColor', 'r', 'LineWidth', 1);
l4 = errorbar(ns, ms10*r^10, vs10*r^20, 'g-s', 'MarkerSize', 4, ...
              'MarkerFaceColor', 'g', 'LineWidth', 1);

h = xlabel('# of samples (n)');
set(h, 'FontSize', 20);
h = ylabel('Absolute Error');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
xlim([ns(1) ns(end)]);
legend([l1 l2 l3 l4], 'd=1', 'd=3', 'd=5', 'd=10');
saveas(f, './figs/conv.fig', 'fig');
saveas(f, './figs/conv.eps', 'epsc');

f = figure();
l1 = plot(ns, sqrt(ns).*ms1*r, 'b-o', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'b', 'LineWidth', 1);
hold on;
l2 = plot(ns, sqrt(ns).*ms3*r^3, 'k-d', 'MarkerSize', 4, ...
          'MarkerFaceColor', 'k', 'LineWidth', 1);
l3 = plot(ns, sqrt(ns).*ms5*r^5, 'r-d', 'MarkerSize', 4, ...
          'MarkerFaceColor', 'r', 'LineWidth', 1);
l4 = plot(ns, sqrt(ns).*ms10*r^10, 'g-s', 'MarkerSize', 4, ...
          'MarkerFaceColor', 'g', 'LineWidth', 1);

h = xlabel('# of samples (n)');
set(h, 'FontSize', 20);
h = ylabel('sqrt(n) * Absolute Error');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
legend([l1 l2 l3 l4], 'd=1', 'd=3', 'd=5', 'd=10', 2);
saveas(f, './figs/conv_rescaled.fig', 'fig');
saveas(f, './figs/conv_rescaled.eps', 'epsc');