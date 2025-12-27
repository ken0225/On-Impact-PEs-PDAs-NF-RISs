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

%% Clean All & Timer
close all;
clear;
tic;
warning off;

%% System Parameters Initialization
disp('Step: Parameter initialization.');
global c_0 f_c lambda_c; 

c_0 = physconst('Lightspeed'); 

f_c = 2.4e9; 

lambda_c = c_0/f_c;

kappa_t =0; kappa_r = kappa_t; 

tau = pi/2; 
mu = 0;
kappa = 2; 

realization = 1e2; 

varrho = 1e3; 

iota = 0; 

sqrt_M = 2e2; 

d_x = lambda_c / 0.5; d_y = d_x; 

M = sqrt_M^2;

r_RIS = 10; 

r_User = 1.5; 

r_AP = 15; 

v_speed = 1; 

Pt_dBm = 20; 

Pt_W = 10^(Pt_dBm/10)/1e3;

sigma_square_dBm = -80; 

sigma_square_W = 10^(sigma_square_dBm/10)/1e3;

x_User = 20; 
y_User = 0+r_User;
z_User = 8; 

x_AP = -x_User;
y_AP = 0+r_AP;
z_AP = z_User;

D_AP = [x_AP, y_AP, z_AP];

D_User = [x_User, y_User, z_User];

%% Step: Compute coordinate and the time delay.
disp('------------Step: Compute coordinate and the time delay.------------');

[centers_RIS]=function_centers_single_RIS(sqrt_M, sqrt_M, d_x, d_y, r_RIS);

x_RIS = 0; 

T = length(x_RIS);

for kk = 1 : T
    
    temp_D_RIS_centers = centers_RIS;
    
    temp_D_RIS_centers(:,1) = centers_RIS(:,1) + x_RIS(kk);
    
    D_RIS_centers(:, :, kk) = temp_D_RIS_centers;
    
    clear temp_D_RIS_centers
    
end

[tau_0, tau_m, tau_AP_RIS, tau_RIS_User, d_0, d_AP_RIS, d_RIS_User]...
    =function_time_delay(D_RIS_centers, D_AP, D_User);

A_0 = 0;

%% Compute dynamic rician factors and optimal RIS phase shift.
disp('------------Step: Compute optimal RIS phase shifts.------------');

[kappa_d, rho_m, kappa_m] = function_Rician_factor(d_0, d_AP_RIS, d_RIS_User, varrho, iota);

K_1 = sqrt(kappa_d ./ (kappa_d+1)); 
K_2 = sqrt(1 ./ (kappa_d+1)); 

K_3 = sqrt(mean(kappa_m, 2) ./ (mean(kappa_m, 2)+1)); 
K_4 = sqrt(1 ./ (mean(kappa_m, 2)+1)); 

G_1 = sqrt(mean(rho_m, 2) ./ (mean(rho_m, 2)+1)); 
G_2 = sqrt(1 ./ (mean(rho_m, 2)+1)); 

zeta_m = ceil((f_c*(tau_m - tau_0)));

phi_m = 2*pi*(f_c*(tau_0 - tau_m)+zeta_m);

phi_m = mod(phi_m + pi, 2*pi) - pi; 

exp_phi_m_matrix = exp(-1i.*phi_m);

exp_phi_m = [];

for ii = 1 : realization
    
    exp_phi_m(:, :, ii) = exp_phi_m_matrix;
    
end

%% Compute projected aperature
disp('------------Step: Compute projected aperature.------------');

pixel_area = d_x * d_y;

projected_angle_RIS_to_AP = [];

for ii = 1 : T
    
    [temp_projected_angle]=function_projected_angle(D_RIS_centers(:, :, ii), D_AP);
    
    projected_angle_RIS_to_AP(:, ii) = temp_projected_angle;
    
end

projected_angle_RIS_to_User = [];

for jj = 1 : T
    
    [temp_projected_angle]=function_projected_angle(D_RIS_centers(:, :, jj), D_User);
    
    projected_angle_RIS_to_User(:, jj) = temp_projected_angle;
      
end

projected_factor_RIS_to_AP = 4*cos(projected_angle_RIS_to_AP);
    
projected_factor_RIS_to_AP(projected_factor_RIS_to_AP<0) = 0;

projected_aperture_RIS_to_AP = pixel_area * projected_factor_RIS_to_AP; 

projected_factor_RIS_to_User = 4*cos(projected_angle_RIS_to_User);
    
projected_factor_RIS_to_User(projected_factor_RIS_to_User<0) = 0;

projected_aperture_RIS_to_User = pixel_area * projected_factor_RIS_to_User;

%% Step: Generate the Rician channel
disp('------------Step: Generate the Rician channel:------------');

fprintf('         Generate the Rician channel. Progress:   0%%');

h_d_NLoS = sqrt(1/2) * (randn(1, realization) + 1i*randn(1, realization));
fprintf('\r         Generate the Rician channel. Progress:  10%%');

g_m_NLoS = sqrt(1/2) * (randn(M, T, realization) + 1i*randn(M, T, realization));
fprintf('\r         Generate the Rician channel. Progress:  30%%');

h_m_NLoS = sqrt(1/2) * (randn(M, T, realization) + 1i*randn(M, T, realization));
fprintf('\r         Generate the Rician channel. Progress:  50%%');

g_m_LoS = exp(-1i*2*pi*f_c*tau_AP_RIS);
fprintf('\r         Generate the Rician channel. Progress:  70%%');

h_m_LoS = exp(-1i*2*pi*f_c*tau_RIS_User);
fprintf('\r         Generate the Rician channel. Progress:  90%%');

fprintf('\r         Generate the Rician channel. Progress:  100%%\r');

%% For PE, we consider uniform and Von Mises RVs.
disp('------------Step: Generate UF and VM RVs:------------');

if ~exist('realization','var') || isempty(realization)
    error('Variable ''realization'' is not defined. Please define it (e.g., realization = 100).');
end

progress_steps = 0:20:100;
next_step_index = 1;
fprintf('         Generate UF noise (delta_m): 0%%');

delta_m = random('Uniform', -tau, tau, T, M, realization);
fprintf('\r         Generate UF noise (bar_delta_m): 20%%');

bar_delta_m = random('Uniform', -tau, tau, T, M, realization);
fprintf('\r         Generate VM noise (gamma_m): 40%%');

gamma_m = function_vmrand(mu, kappa, T, M, realization);
fprintf('\r         Generate VM noise (bar_gamma_m): 60%%');

bar_gamma_m = function_vmrand(mu, kappa, T, M, realization);
fprintf('\r         Generate PE: 80%%');
fprintf('\r         Generate PE: 100%%\n');

%% SE Calculation
disp('------------Step: Analytical Results SE calculation.------------');

a=1; b=0.2; c=0.43*pi;

Gamma_SE = (((1-b)/2 .* sin(phi_m-c) .* (sin(tau)/tau) .* (besseli(1,kappa)./besseli(0,kappa)) + (1+b)/2) .* (besseli(1,kappa)./besseli(0,kappa)) .* (sin(tau)/tau)).^2;

h_abs_square_SE = (d_x^2*d_y^2*D_AP(3)*D_User(3))/(pi^2) .* sum(sqrt(Gamma_SE)./(d_AP_RIS.*d_RIS_User).^1.5, 2).^2;

SNR = Pt_W*h_abs_square_SE/sigma_square_W;

SE_AN = log2(1+SNR);

%% Step: Plot the phase distribution
figure(1); 
edges = linspace(-pi, pi, 100); 

phi_m_flattened = phi_m(:); 

pdf_values = histogram(phi_m_flattened, edges, 'Normalization', 'pdf');

bin_centers = edges(1:end-1) + diff(edges)/2;
pdf_data = pdf_values.Values;

area(bin_centers, pdf_data, 'FaceColor', [169, 209, 142] / 255, 'EdgeColor', 'none'); 

hold on; 

plot(bin_centers, pdf_data, 'LineWidth', 3, 'Color', [244, 177, 131] / 255); 
xlim([-pi, pi]);

xlabel('$\phi_m\in[-\pi, \pi]$', 'Interpreter', 'LaTex');
ylabel('Normalized PDF', 'Interpreter', 'LaTeX');

xticks(-pi:pi/2:pi); 
xticklabels({'$-\pi$', '$-{\pi}/{2}$', '$0$', '${\pi}/{2}$', '$\pi$'}); 
ax = gca; 
ax.TickLabelInterpreter = 'latex'; 

set(gca,'FontSize',16, 'LineWidth', 1.2);

grid on;

hold off;

%% Step: Plot the Gamma_SE distribution
figure(2); 

Gamma_SE_flattened = Gamma_SE(:); 

qqq = 0.10;  
num_greater = sum(Gamma_SE_flattened > qqq);
N = numel(Gamma_SE_flattened);
fprintf('Gamma_SE > %.4f percentage: %.2f%% (%d / %d)\n', qqq, num_greater / N * 100, num_greater, N);
fprintf('Gamma_SE < %.4f percentage: %.2f%% (%d / %d)\n', qqq, 100-num_greater / N * 100, N-num_greater, N);

rrr = 0.1;   
www = 0.12;   
if rrr > www
    temp = rrr;
    rrr = www;
    www = temp;
end
num_in_range = sum(Gamma_SE_flattened >= rrr & Gamma_SE_flattened <= www);
N = numel(Gamma_SE_flattened);
percentage = num_in_range / N * 100;
fprintf('Gamma_SE in [%.4f, %.4f] samples: %d / %d, percentage: %.2f%%\n', ...
    rrr, www, num_in_range, N, percentage);

edges = linspace(min(Gamma_SE_flattened), max(Gamma_SE_flattened), 100);

pdf_values = histogram(Gamma_SE_flattened, edges, 'Normalization', 'pdf');

bin_centers = edges(1:end-1) + diff(edges)/2;
pdf_data = pdf_values.Values;

pdf_data_normalized = pdf_data / max(pdf_data);

area(bin_centers, pdf_data_normalized, 'FaceColor', [169, 209, 142] / 255, 'EdgeColor', 'none'); 

hold on; 

plot(bin_centers, pdf_data_normalized, 'LineWidth', 3, 'Color', [244, 177, 131] / 255); 

xlabel('$\Gamma$', 'Interpreter', 'LaTex');
ylabel('Normalized PDF', 'Interpreter', 'LaTeX');
xlim([min(Gamma_SE_flattened)-0.001, max(Gamma_SE_flattened)+0.001]);
ylim([0, 1]); 

set(gca,'FontSize',16, 'LineWidth', 1.2);
ax = gca;
ax.TickLabelInterpreter = 'latex';   
set(gca,'FontSize',16, 'LineWidth', 1.2);

grid on;

hold off;

%% Step: UF pi/2, pi/4, pi/8 Gamma_SE PDF Comparison
tau_sets = [pi/2, pi/4, pi/8]; 
num_angles = length(tau_sets);
Gamma_SE_results = cell(num_angles, 1);

for index = 1:num_angles
    tau = tau_sets(index); 
    Gamma_SE = (((1-b)/2 .* sin(phi_m-c) .* (sin(tau)/tau) .* (besseli(1,kappa)./besseli(0,kappa)) + (1+b)/2) .* (besseli(1,kappa)./besseli(0,kappa)) .* (sin(tau)/tau)).^2;
    Gamma_SE_results{index} = Gamma_SE;
end

figure(3); 
hold on; 

colors = {[244, 177, 131] / 255,  [169, 209, 142] / 255, [100, 149, 237] ./ 255};
angle_labels = {'$\tau = {\pi}/{2}$', '$\tau = {\pi}/{4}$', '$\tau = {\pi}/{8}$'};

num_bins = 100;
bin_centers_all = zeros(num_angles, num_bins);
pdf_data_normalized_all = zeros(num_angles, num_bins);

for index = 1:num_angles
    Gamma_SE_flattened = Gamma_SE_results{index}(:);
    edges = linspace(min(Gamma_SE_flattened), max(Gamma_SE_flattened), num_bins + 1);
    [pdf_data, ~] = histcounts(Gamma_SE_flattened, edges, 'Normalization', 'pdf');
    bin_centers = edges(1:end-1) + diff(edges)/2;
    
    if max(pdf_data) > 0
        pdf_data_normalized = pdf_data / max(pdf_data);
    else
        pdf_data_normalized = pdf_data;
    end
    
    bin_centers_all(index, :) = bin_centers;
    pdf_data_normalized_all(index, :) = pdf_data_normalized;
end

for index = 1:num_angles
    plot(bin_centers_all(index, :), pdf_data_normalized_all(index, :), 'LineWidth', 3, 'Color', colors{index});
end

xlabel('$\Gamma$', 'Interpreter', 'LaTeX');
ylabel('Normalized PDF Envelope', 'Interpreter', 'LaTeX');

ylim([0, 1]); 

legend(angle_labels, 'Interpreter', 'LaTeX', 'Location', 'Best');

ax = gca;
ax.TickLabelInterpreter = 'latex';   
set(gca,'FontSize',16, 'LineWidth', 1.2);

grid on; box on;

hold off;

%% Timer ends
toc;