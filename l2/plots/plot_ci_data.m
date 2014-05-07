
load('./c_data_1_2_3.mat');

f = figure();
l1 = plot(ns, ps1, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
l2 = plot(ns, ps2,  'r-d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
l3 = plot(ns, ps3, 'g-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');

h = xlabel('# of samples (n)');
set(h, 'FontSize', 20);
h = ylabel('Success Probability');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
legend([l1 l2 l3], 'd=1', 'd=2', 'd=3', 4);
saveas(f, './figs/ci_trap.fig', 'fig');
saveas(f, './figs/ci_trap.eps', 'epsc');
