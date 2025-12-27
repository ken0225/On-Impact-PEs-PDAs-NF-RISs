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

close all;
clear;

% Range of frequency values
fc = 2.4;                        % [GHz]
% Compute angular frequency
omega = 2*pi*fc*1e9;             % [rad/s]

% Varying capacitance value
Cn = linspace(0.47, 2.35, 1000)*1e-12;   % [F]
Rn = [0.5, 1, 1.5, 2, 2.5];              % [Ohm]

% Prepare to store results
vn = zeros(length(Rn), length(Cn));
for n = 1:length(Rn)
    for m = 1:length(Cn)
        vn(n, m) = function_refcoefficient(omega, Cn(m), Rn(n));
    end
end

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

%% Proposed phase shift model
b = 0.2;
a = 1.6;
c = 0.43*pi;

% Sparse points to obtain uniformly distributed dots
num_points_beta = 30;
theta_mn = linspace(-pi, pi, num_points_beta);
beta_mn  = zeros(1, length(theta_mn));
for ll = 1:length(theta_mn)
    beta_mn(1, ll) = function_beta_mn(theta_mn(ll), a, b, c);
end

%% Plot
figure;
hold on; box on; grid on;

yyaxis left;
plot(angle(vn(1,:)), abs(vn(1,:)), '-', 'LineWidth', 3.5, 'Color', [169,209,142]./255);
plot(angle(vn(2,:)), abs(vn(2,:)), '-', 'LineWidth', 3.5, 'Color', [169,209,142]./255);
plot(angle(vn(3,:)), abs(vn(3,:)), '-', 'LineWidth', 3.5, 'Color', [169,209,142]./255);
plot(angle(vn(4,:)), abs(vn(4,:)), '-', 'LineWidth', 3.5, 'Color', [169,209,142]./255);
pL = plot(angle(vn(5,:)), abs(vn(5,:)), '-', 'LineWidth', 3.5, 'Color', [169,209,142]./255);
ylabel('$\varsigma_m$','Interpreter','latex');
set(gca, 'YColor', 'k');   % Set left y-axis color to black
ylim([b-0.05, 1+0.03]);    

yyaxis right;
pR = plot(theta_mn, beta_mn, 'o', 'LineWidth', 2, 'MarkerSize', 12, ...
          'Color', [244,177,131]./255);
ylabel('$\beta(\phi_{m})$','Interpreter','latex');
set(gca, 'YColor', 'k');   % Set right y-axis color to black
ylim([b-0.05, 1+0.03]);

xlabel('$\phi_{m}\in[-{\pi}/{2}-c,\, {\pi}/{2}+c]$','Interpreter','latex');

% Set x-axis range
xlim([-pi-0.03, pi+0.03]);
title('$a=1.6,\,b=0.2, \,\mathrm{and}\ c=0.43\pi$', 'Interpreter', 'latex');

% Display ticks as fractions of pi
xticks([-pi, -pi/2, 0, pi/2, pi]);
xticklabels({ ...
    '$-\pi$', '$-\tfrac{\pi}{2}$', '$0$', '$\tfrac{\pi}{2}$', '$\pi$' ...
});

% Adjust view and label angles
view(2);                 % 2D top view
xtickangle(0);           % X-axis tick angle
ytickangle(0);           % Y-axis tick angle
set(gca, 'TickDir', 'out');  % Ticks point outwards

legend([pL(1), pR(1)], ...
    'Reflection Coefficient $\varsigma_m$', ...
    'Approximated PDA $\beta(\phi_m)$', ...
    'Interpreter', 'latex', 'Location', 'northeast');

text(0.6, 0.185, '$R_m=0.5$ to $2.5 \Omega$', ...
     'Interpreter', 'latex', 'Color', 'Black', 'FontSize', 16);

set(gca, 'FontSize', 16, 'LineWidth', 1.2);