clear
close all
load converge_revise15.mat

y = var_change(1:end, 1);
x = 1:length(y);

figure
plot(x, y, 'b-', 'linewidth', 2);
set(gca, 'fontsize', 40, 'linewidth', 4);
xlabel(' ', 'FontSize', 40);
ylabel('Change in design variables', 'FontSize', 40);
hold on
y1 = var_change(1:end, 2);
hold on
plot(x, y1, 'r-', 'linewidth', 2);
% axis([0, 500 -1 30])

legend('Change in \phi', 'Change in \mu', 'fontsize', 30);