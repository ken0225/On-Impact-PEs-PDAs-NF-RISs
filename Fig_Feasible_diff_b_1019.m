%% Program Information
% Author: Ke (Ken) Wang, Macao Polytechnic University
% Email: ke.wang@mpu.edu.mo; kewang0225@gmail.com
% Personal Page: https://kewang.fun
%
% License: GNU General Public License v2.0 (GPLv2)
%
% Citation:
% If you use this code in research that leads to publications, please cite:
% K. Wang, C.-T. Lam, B. K. Ng, and Y. Liu, "On the Impact of Phase Errors in
% Phase-Dependent Amplitudes of Near-Field RISs," IEEE Transactions on
% Vehicular Technology, doi: 10.1109/TVT.2025.3647594
%
% Link: https://ieeexplore.ieee.org/document/11314716

%% Practical Phase Shift Model Visualization (Varying b values)
% Generates polar plots for different minimum amplitude b values

clear; clc; close all;

% Define angle range
theta = 0:0.01:2*pi;

% Parameter settings
a = 1;
c = 0.43*pi;
b_values = [0.2, 0.4, 0.6, 0.8];
% Colors based on previous reference style
colors = {[169, 209, 142]./255, [244, 177, 131]./255, [231, 98, 84]./255, [128, 116, 200]./255}; 

% Create figure
figure;
hold on; grid on;

% Plot ideal phase shift model (Unit Circle)
ideal_circle = exp(1j * theta);
plot(real(ideal_circle), imag(ideal_circle), 'k:', 'LineWidth', 3, 'DisplayName', 'Ideal Unit Circle');

% Loop through different 'b' values to plot curves
for i = 1:length(b_values)
    b = b_values(i);
    alpha = (1 - b) * ((sin(theta - c) + 1)/2).^a + b;
    z = alpha .* exp(1j * theta);
    
    % Calculate the area of the convex hull
    hull_indices = convhull(real(z), imag(z));
    area_z = polyarea(real(z(hull_indices)), imag(z(hull_indices)));
    
    % Plot curve with specific RGB color and display area in legend
    plot(real(z), imag(z), 'Color', colors{i}, 'LineWidth', 3, ...
        'DisplayName', sprintf('$b=%.1f$ (Area = %.3f)', b, area_z));
end

% Plot origin point
plot(0, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');

% --- Uniform style graphical attribute settings ---

% Set axes limits and labels
axis normal;
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
xlabel('$\mathsf{Re}$', 'Interpreter', 'latex');
ylabel('$\mathsf{Im}$', 'Interpreter', 'latex');

% Add title showing fixed parameters
title(sprintf('$a=%.0f, c=%.2f\\pi$', a, c/pi), 'Interpreter', 'latex');

% Configure legend and axis font properties
legend('Location', 'SouthEast', 'Interpreter', 'latex');
ax = gca;
set(ax, 'XTick', -1.5:0.5:1.5);
set(ax, 'YTick', -1.5:0.5:1.5);
ax.TickLabelInterpreter = 'latex';   % Applies to both axis tick labels
set(ax, 'FontSize', 16, 'LineWidth', 1.2);
box on;
hold off;