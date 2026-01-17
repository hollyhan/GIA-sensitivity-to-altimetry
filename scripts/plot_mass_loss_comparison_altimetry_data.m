% plot_mass_loss_comparison.m
% Script to plot total mass loss comparison across different altimetry datasets

% Data from the chart
datasets = {'JPL-GEMB', 'JPL-GSFC', 'DTU2016', 'DTU2025', 'Buffalo-GEMB', 'Buffalo-GSFC', 'Buffalo-IMAU'};
years = {'(1993-2024)', '(1993-2024)', '(2002-2022)', '(2002-2022)', '(1994-2020)', '(1994-2020)', '(1994-2020)'};
mass_loss_gt = [5177, 5899, 4165, 4191, 4899, 5199, 4879];

% Create figure
figure('Position', [100, 100, 1200, 600]);

% Create bar plot
b = bar(mass_loss_gt, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'black', 'LineWidth', 1.5);

% Customize the plot
set(gca, 'XTickLabel', datasets, 'FontSize', 16);
ylabel('Total Mass Loss (Gt)', 'FontSize', 16, 'FontWeight', 'bold');
title('Total Mass Loss Comparison Across Altimetry Datasets', 'FontSize', 16, 'FontWeight', 'bold');

% Add value labels on top of each bar
for i = 1:length(mass_loss_gt)
    text(i, mass_loss_gt(i) + 100, num2str(mass_loss_gt(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
end

% Add year range labels below x-axis
ax = gca;
ax.XTickLabel = {};
for i = 1:length(datasets)
    text(i, -300, [datasets{i}, newline, years{i}], ...
         'HorizontalAlignment', 'center', 'FontSize', 16, 'Rotation', 0);
end

% Adjust y-axis limits to accommodate labels
ylim([0, max(mass_loss_gt) + 500]);

% Add grid
grid on;
grid minor;

% Customize appearance
set(gca, 'FontSize', 16, 'LineWidth', 1.2);
box on;

% Add legend for color coding (optional)
legend_entries = {'Mass Loss (Gt)'};
legend(legend_entries, 'Location', 'northeast', 'FontSize', 12);

% Add some statistics
mean_loss = mean(mass_loss_gt);
std_loss = std(mass_loss_gt);
text(0.02, 0.98, sprintf('Mean: %.0f ± %.0f Gt', mean_loss, std_loss), ...
     'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', ...
     'VerticalAlignment', 'top', 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Save the figure
saveas(gcf, 'mass_loss_comparison.png', 'png');
saveas(gcf, 'mass_loss_comparison.pdf', 'pdf');

fprintf('Mass loss comparison plot saved as mass_loss_comparison.png and mass_loss_comparison.pdf\n');
fprintf('Mean mass loss: %.0f ± %.0f Gt\n', mean_loss, std_loss); 