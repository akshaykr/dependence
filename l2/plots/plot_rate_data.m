load('./rate_data.mat');

r = 2*sqrt(pi); %% l_2^2 of 1-d gaussian density.

f = figure();
l1 = errorbar(ns, ms1*r, vs1*r^2, 'b-o', 'MarkerSize', 4, ...
              'MarkerFaceColor', 'b', 'LineWidth', 1);
hold on;
l2 = errorbar(ns, ms2*r^2, vs2*r^4, 'r-d', 'MarkerSize', 4, ...
              'MarkerFaceColor', 'r', 'LineWidth', 1);
l3 = errorbar(ns, ms3*r^3, vs3*r^6, 'g-s', 'MarkerSize', 4, ...
              'MarkerFaceColor', 'g', 'LineWidth', 1);

h = xlabel('# of samples (n)');
set(h, 'FontSize', 20);
h = ylabel('Absolute Error');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
legend([l1 l2 l3], 'd=1', 'd=2', 'd=3');
saveas(f, './figs/conv.fig', 'fig');
saveas(f, './figs/conv.eps', 'epsc');

f = figure();
l1 = plot(ns, sqrt(ns).*ms1*r, 'b-o', 'MarkerSize', 4, 'MarkerFaceColor', ...
          'b', 'LineWidth', 1);
hold on;
l2 = plot(ns, sqrt(ns).*ms2*r^2, 'r-d', 'MarkerSize', 4, ...
          'MarkerFaceColor', 'r', 'LineWidth', 1);
l3 = plot(ns, sqrt(ns).*ms3*r^3, 'g-s', 'MarkerSize', 4, ...
          'MarkerFaceColor', 'g', 'LineWidth', 1);

h = xlabel('# of samples (n)');
set(h, 'FontSize', 20);
h = ylabel('sqrt(n) * Absolute Error');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
legend([l1 l2 l3], 'd=1', 'd=2', 'd=3');
saveas(f, './figs/conv_rescaled.fig', 'fig');
saveas(f, './figs/conv_rescaled.eps', 'epsc');