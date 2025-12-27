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

%% Initialization
close all; clear; tic; warning off;

%% System Parameters
disp('Parameter initialization.');
global c_0 f_c lambda_c;
c_0 = physconst('Lightspeed');
f_c = 2.4e9;
lambda_c = c_0 / f_c;

% RIS Configuration
sqrt_M = 2e2;
d_x = lambda_c / 2;
d_y = d_x;
M = sqrt_M^2;

% Coordinates
r_RIS = 10;
r_User = 1.5;
r_AP = 15;
x_User = 20;
y_User = r_User;
z_User = 8;
x_AP = -x_User;
y_AP = r_AP;
z_AP = z_User;

D_AP = [x_AP, y_AP, z_AP];
D_User = [x_User, y_User, z_User];

%% Coordinate and Time Delay Computation
disp('Compute coordinate and time delay.');

% RIS Element Centers
[D_RIS_centers] = function_centers_single_RIS(sqrt_M, sqrt_M, d_x, d_y, r_RIS);

%% Common Gamma Constants
disp('Compute Gamma constants.');

b = 0.2;
c = 0.43 * pi;

%% phi_m Calculation (Only Once)
disp('Compute phi_m values.');

% phi_m Preallocation
phi_m_values = zeros(sqrt_M, sqrt_M);

% phi_m Computation
for ii = 1:sqrt_M
    for jj = 1:sqrt_M
        % Coordinates
        x_coord = D_RIS_centers(ii, 1);
        y_coord = D_RIS_centers(ii, 2);
        D_RIS_temp = [x_coord, y_coord, r_RIS];

        % Time Delays
        [tau_0_temp, tau_m_temp] = function_time_delay(D_RIS_temp, D_AP, D_User);

        % Phase Shift
        zeta_m_temp = ceil(f_c * (tau_m_temp - tau_0_temp));
        phi_m_values(ii, jj) = mod(2 * pi * (f_c * (tau_0_temp - tau_m_temp) + zeta_m_temp) + pi, 2 * pi) - pi;
    end
end

% Find min and max phi_m values
min_phi_m = min(phi_m_values(:));
max_phi_m = max(phi_m_values(:));

%% Gamma Calculation for Different kappa Values at tau = pi/8
tau = pi/8; % Fixed tau value
kappa_values = [2, 5, 8];

% Preallocate to store Gamma values
Gamma_all = cell(size(kappa_values));

% Calculate Gamma values for each kappa
for k = 1:length(kappa_values)
    kappa = kappa_values(k); % Current kappa value

    % Gamma Preallocation
    Gamma_values = zeros(sqrt_M, sqrt_M);

    % Gamma Computation
    for ii = 1:sqrt_M
        for jj = 1:sqrt_M
            % Coordinates
            x_coord = D_RIS_centers(ii, 1);
            y_coord = D_RIS_centers(ii, 2);
            D_RIS_temp = [x_coord, y_coord, r_RIS];

            % Time Delays
            [tau_0_temp, tau_m_temp] = function_time_delay(D_RIS_temp, D_AP, D_User);

            % Gamma Value
            Gamma_values(ii, jj) = (((1 - b) / 2 .* sin(phi_m_values(ii, jj) - c) .* (sin(tau) / tau) .* ...
                (besseli(1, kappa) ./ besseli(0, kappa)) + (1 + b) / 2) .* ...
                (besseli(1, kappa) ./ besseli(0, kappa)) .* (sin(tau) / tau)).^2;
        end
    end

    % Store Gamma values
    Gamma_all{k} = Gamma_values;
end

% Find min and max Gamma values across all kappa values
min_Gamma = inf;
max_Gamma = -inf;

for k = 1:length(kappa_values)
    min_Gamma = min(min_Gamma, min(Gamma_all{k}(:)));
    max_Gamma = max(max_Gamma, max(Gamma_all{k}(:)));
end

% %% Plotting - phi_m Heatmap (Only Once)
% figure(1);  % Create a new figure for phi_m
% 
% imagesc(phi_m_values);
% set(gca,'CLim',[-pi, pi]);
% 
% hc = colorbar;
% set(hc,'TickLabelInterpreter','latex',...
%        'Ticks',[-pi -pi/2 0 pi/2 pi],...
%        'TickLabels',{'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'});
% colormap(parula);
% 
% % Axis Labels and Title with LaTeX
% axis equal; % Ensure the plot is square
% axis([0, 200, 0, 200]); % Setting axis limits after axis equal
% xlabel('$\sqrt{M}\in[0,200]$', 'Interpreter', 'latex');
% ylabel('$\sqrt{M}\in[0,200]$', 'Interpreter', 'latex');
% title('$\phi_m\in[-\pi, \pi]$', 'Interpreter', 'latex');
% 
% % Font Settings
% set(gca, 'FontSize', 16, 'LineWidth', 1.2);
% grid on;
% 
% %% Plotting - Separate Figures for kappa
% 
% % --- Plot for kappa = 2 ---
% figure(5);
% imagesc(Gamma_all{1});
% clim([0, 1]);
% colorbar;
% colormap(parula);
% axis equal;
% axis([0, 200, 0, 200]);
% xlabel('$\sqrt{M}\in[0,200]$', 'Interpreter', 'latex');
% ylabel('$\sqrt{M}\in[0,200]$', 'Interpreter', 'latex');
% title('$\tau = \pi/8$, $\kappa = 2$', 'Interpreter', 'latex');
% set(gca, 'FontSize', 16, 'LineWidth', 1.2');
% grid on;
% 
% % --- Plot for kappa = 5 ---
% figure(6);
% imagesc(Gamma_all{2});
% clim([0, 1]);
% colorbar;
% colormap(parula);
% axis equal;
% axis([0, 200, 0, 200]);
% xlabel('$\sqrt{M}\in[0,200]$', 'Interpreter', 'latex');
% ylabel('$\sqrt{M}\in[0,200]$', 'Interpreter', 'latex');
% title('$\tau = \pi/8$, $\kappa = 5$', 'Interpreter', 'latex');
% set(gca, 'FontSize', 16, 'LineWidth', 1.2);
% grid on;
% 
% % --- Plot for kappa = 8 ---
% figure(7);
% imagesc(Gamma_all{3});
% clim([0, 1]);
% colorbar;
% colormap(parula);
% axis equal;
% axis([0, 200, 0, 200]);
% xlabel('$\sqrt{M}\in[0,200]$', 'Interpreter', 'latex');
% ylabel('$\sqrt{M}\in[0,200]$', 'Interpreter', 'latex');
% title('$\tau = \pi/8$, $\kappa = 8$', 'Interpreter', 'latex');
% set(gca, 'FontSize', 16, 'LineWidth', 1.2);
% grid on;

%% Plotting - phi_m Heatmap (Only Once)
figure(1);  % Create a new figure for phi_m

imagesc(phi_m_values);
set(gca,'CLim',[-pi, pi]);

hc = colorbar;
set(hc,'TickLabelInterpreter','latex',...
       'Ticks',[-pi -pi/2 0 pi/2 pi],...
       'TickLabels',{'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'});
colormap(parula);

% Axis Labels and Title with LaTeX
axis equal; % Ensure the plot is square
axis([0, 200, 0, 200]); % Setting axis limits after axis equal
xlabel('$\sqrt{M}\in[0,200]$','Interpreter','latex');
ylabel('$\sqrt{M}\in[0,200]$','Interpreter','latex');
title('$\phi_m\in[-\pi,\pi]$','Interpreter','latex');

% Font Settings
set(gca,'FontSize',16,'LineWidth',1.2);
grid on;

%% Plotting - Separate Figures for kappa

% --- Plot for kappa = 2 ---
figure(5);
imagesc(Gamma_all{1});
clim([0, 1]);
hc2 = colorbar;
set(hc2,'TickLabelInterpreter','latex');   % 仅启用 LaTeX，不改你的数值范围
colormap(parula);
axis equal;
axis([0, 200, 0, 200]);
xlabel('$\sqrt{M}\in[0,200]$','Interpreter','latex');
ylabel('$\sqrt{M}\in[0,200]$','Interpreter','latex');
title('$\tau=\pi/8,\ \kappa=2$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');  % 主图坐标轴刻度用 LaTeX
set(gca,'FontSize',16,'LineWidth',1.2);
grid on;

% --- Plot for kappa = 5 ---
figure(6);
imagesc(Gamma_all{2});
clim([0, 1]);
hc3 = colorbar;
set(hc3,'TickLabelInterpreter','latex');
colormap(parula);
axis equal;
axis([0, 200, 0, 200]);
xlabel('$\sqrt{M}\in[0,200]$','Interpreter','latex');
ylabel('$\sqrt{M}\in[0,200]$','Interpreter','latex');
title('$\tau=\pi/8,\ \kappa=5$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');  % 主图坐标轴刻度用 LaTeX
set(gca,'FontSize',16,'LineWidth',1.2);
grid on;

% --- Plot for kappa = 8 ---
figure(7);
imagesc(Gamma_all{3});
clim([0, 1]);
hc4 = colorbar;
set(hc4,'TickLabelInterpreter','latex');
colormap(parula);
axis equal;
axis([0, 200, 0, 200]);
xlabel('$\sqrt{M}\in[0,200]$','Interpreter','latex');
ylabel('$\sqrt{M}\in[0,200]$','Interpreter','latex');
title('$\tau=\pi/8,\ \kappa=8$','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');  % 主图坐标轴刻度用 LaTeX
set(gca,'FontSize',16,'LineWidth',1.2);
grid on;


%% End
toc;
