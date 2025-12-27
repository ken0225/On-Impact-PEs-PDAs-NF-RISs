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

%% Parameter Sweep: Different b and a values
close all; clear;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

% Angular sampling (smooth enough for plotting)
numPts = 1000;                          
theta_mn = linspace(-pi, pi, numPts);

%% different b
b_all = [0.2, 0.4, 0.6, 0.8, 1.0];
a_fix = 1.6;
c = 0.43*pi;

beta_mn_all_b = zeros(length(b_all), numPts);
for kk = 1:length(b_all)
    % Vectorized calculation using arrayfun
    beta_mn_all_b(kk, :) = arrayfun(@(t) function_beta_mn(t, a_fix, b_all(kk), c), theta_mn);
end

%% different a
a_all = [1.0, 1.4, 1.6, 2.0, 2.4];
b_fix = 0.2;

beta_mn_all_a = zeros(length(a_all), numPts);
for jj = 1:length(a_all)
    beta_mn_all_a(jj, :) = arrayfun(@(t) function_beta_mn(t, a_all(jj), b_fix, c), theta_mn);
end

%% Plotting
figure; hold on; box on; grid on;

% Different a (Orange, dash-dot line)
colA = [244, 177, 131]./255;
p6 = plot(theta_mn, beta_mn_all_a(1, :), '-.', 'LineWidth', 3.5, 'Color', colA);
p7 = plot(theta_mn, beta_mn_all_a(2, :), '-.', 'LineWidth', 3.5, 'Color', colA);
p8 = plot(theta_mn, beta_mn_all_a(3, :), '-.', 'LineWidth', 3.5, 'Color', colA);
p9 = plot(theta_mn, beta_mn_all_a(4, :), '-.', 'LineWidth', 3.5, 'Color', colA);
p10= plot(theta_mn, beta_mn_all_a(5, :), '-.', 'LineWidth', 3.5, 'Color', colA);

% Different b (Green, solid line)
colB = [169, 209, 142]./255;
p1 = plot(theta_mn, beta_mn_all_b(1, :), '-', 'LineWidth', 3.5, 'Color', colB);
p2 = plot(theta_mn, beta_mn_all_b(2, :), '-', 'LineWidth', 3.5, 'Color', colB);
p3 = plot(theta_mn, beta_mn_all_b(3, :), '-', 'LineWidth', 3.5, 'Color', colB);
p4 = plot(theta_mn, beta_mn_all_b(4, :), '-', 'LineWidth', 3.5, 'Color', colB);
p5 = plot(theta_mn, beta_mn_all_b(5, :), '-', 'LineWidth', 3.5, 'Color', colB);

% Text annotations
text(-1.2, 0.5, '$a = 1$ to $2.4$', 'Interpreter', 'latex', 'Color', 'Black', 'FontSize', 16);
text( 1.4, 0.3, '$b = 0.2$ to $1$',  'Interpreter', 'latex', 'Color', 'Black', 'FontSize', 16);

% Axis and labels
xlabel('$\phi_{m}\in[-{\pi}/{2}-c,\, {\pi}/{2}+c]$','Interpreter','latex');
ylabel('$\beta(\phi_{m})$', 'Interpreter', 'latex');
title('$c=0.43\pi$', 'Interpreter', 'latex');

% Set axis limits
axis([-pi-0.03, pi+0.03, b_fix-0.05, 1+0.05]);

% X-axis with pi-fraction ticks
xticks([-pi, -3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi]);
xticklabels({ ...
    '$-\pi$', '$-\tfrac{3\pi}{4}$', '$-\tfrac{\pi}{2}$', '$-\tfrac{\pi}{4}$', ...
    '$0$', '$\tfrac{\pi}{4}$', '$\tfrac{\pi}{2}$', '$\tfrac{3\pi}{4}$', '$\pi$' ...
});

% View configuration and orientation correction
view(2);
camproj('orthographic');
xtickangle(0);
ytickangle(0);

legend([p6(1), p1(1)], ...
    'Different $a$, $b\equiv 0.2$', ...
    'Different $b$, $a\equiv 1.6$', ...
    'Interpreter', 'latex', 'Location', 'northeast');

set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'TickDir', 'out');