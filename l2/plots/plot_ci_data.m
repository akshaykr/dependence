
%% load('./c_data_1_2_3.mat');
load('./ci_data_1.mat');
load('./ci_data_2.mat');
load('./ci_data_3.mat');
load('./ci_data_5.mat');

f = figure();
l1 = plot(ns, ps1, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
l2 = plot(ns, ps2,  'r-d', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
l3 = plot(ns, ps3, 'g-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
l4 = plot(ns, ps5, 'k-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'k');

h = xlabel('# of samples (n)');
set(h, 'FontSize', 20);
h = ylabel('Success Probability');
set(h, 'FontSize', 20);
set(gca, 'FontSize', 20);
%% legend([l1 l2 l3], 'd=1', 'd=2', 'd=3', 4);
legend([l1 l2 l3, l4], 'd=1', 'd=2', 'd=3', 'd=5', 4);
saveas(f, './figs/ci_trap_large.fig', 'fig');
saveas(f, './figs/ci_trap_large.eps', 'epsc');
