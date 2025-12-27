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

%% Parameter Sweep: Different c values
close all; clear;

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

% Angular sampling (smooth enough for plotting)
numPts = 1000;
theta_mn = linspace(-pi, pi, numPts);

%% Fixed parameters
a_fix = 1.6;
b_fix = 0.2;
c_all = [0.5*pi, pi, 1.5*pi, 2*pi];

%% Sweep different c
beta_mn_all_c = zeros(length(c_all), numPts);
for kk = 1:length(c_all)
    % Vectorized calculation for each c value
    beta_mn_all_c(kk, :) = arrayfun(@(t) function_beta_mn(t, a_fix, b_fix, c_all(kk)), theta_mn);
end

%% Plotting
figure; hold on; box on; grid on;

% Different c (using distinct layered colors)
colC = lines(length(c_all));  % Automatically generate distinguishable colors
p = gobjects(length(c_all),1);
for kk = 1:length(c_all)
    p(kk) = plot(theta_mn, beta_mn_all_c(kk, :), '-', 'LineWidth', 3.5, 'Color', colC(kk,:));
end

% Axis and labels
xlabel('$\phi_{m}\in[-{\pi}/{2}-c,\, {\pi}/{2}+c]$', 'Interpreter', 'latex');
ylabel('$\beta(\phi_{m})$', 'Interpreter', 'latex');
title('$a=1.6,\,b=0.2$', 'Interpreter', 'latex');

% Set axis limits
axis([-pi-0.03, pi+0.03, b_fix-0.05, 1+0.05]);

% X-axis with pi-fraction ticks
xticks([-pi, -pi/2, 0, pi/2, pi]);
xticklabels({ ...
    '$-\pi$', '$-\tfrac{\pi}{2}$', '$0$', '$\tfrac{\pi}{2}$', '$\pi$' ...
});

% Visual corrections
view(2);
camproj('orthographic');
xtickangle(0);
ytickangle(0);

% Legend
legend(p, ...
    {'$c=0.5\pi$', '$c=\pi$', '$c=1.5\pi$', '$c=2\pi$'}, ...
    'Interpreter', 'latex', 'Location', 'northeast');

set(gca, 'FontSize', 16, 'LineWidth', 1.2, 'TickDir', 'out');