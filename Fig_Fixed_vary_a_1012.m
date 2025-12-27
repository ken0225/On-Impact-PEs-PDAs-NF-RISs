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
kappa_t = 0; kappa_r = kappa_t;
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

[tau_0, tau_m, tau_AP_RIS, tau_RIS_User, d_0, d_AP_RIS, d_RIS_User] = function_time_delay(D_RIS_centers, D_AP, D_User);
A_0 = 0;

%% Parameter sweep for 'a'
a_values = 1:5; 
num_a = length(a_values);

SE_L_MC_values = zeros(1, num_a);
SE_U_MC_values = zeros(1, num_a);
SE_MC_values = zeros(1, num_a);

%% Main Loop for different 'a' values
for a_idx = 1:num_a
    a = a_values(a_idx);
    disp(['------------ Current a value: ', num2str(a), ' ------------']);

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

    %% Step: Compute projected aperature
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

    fprintf('         Generate the Rician channel. Progress:   0%%\n');

    h_d_NLoS = sqrt(1/2) * (randn(1, realization) + 1i*randn(1, realization));
    fprintf('         Generate the Rician channel. Progress:  10%%\n');

    g_m_NLoS = sqrt(1/2) * (randn(M, T, realization) + 1i*randn(M, T, realization));
    fprintf('         Generate the Rician channel. Progress:  30%%\n');

    h_m_NLoS = sqrt(1/2) * (randn(M, T, realization) + 1i*randn(M, T, realization));
    fprintf('         Generate the Rician channel. Progress:  50%%\n');

    g_m_LoS = exp(-1i*2*pi*f_c*tau_AP_RIS);
    fprintf('         Generate the Rician channel. Progress:  70%%\n');

    h_m_LoS = exp(-1i*2*pi*f_c*tau_RIS_User);
    fprintf('         Generate the Rician channel. Progress:  90%%\n');

    fprintf('         Generate the Rician channel. Progress: 100%%\n');

    %% For PE, we consider uniform and Von Mises RVs.
    disp('------------Step: Generate UF and VM RVs:------------');
    x = pi/2;
    mu = 0;
    kappa = 2;

    if ~exist('realization','var') || isempty(realization)
        error('Variable ''realization'' is not defined. Please define it (e.g., realization = 100).');
    end

    progress_steps = 0:20:100;
    next_step_index = 1;
    fprintf('         Generate UF noise (delta_m): 0%%\n');

    delta_m = random('Uniform', -x, x, T, M, realization);
    fprintf('         Generate UF noise (bar_delta_m): 20%%\n');

    bar_delta_m = random('Uniform', -x, x, T, M, realization);
    fprintf('         Generate VM noise (gamma_m): 40%%\n');

    gamma_m = function_vmrand(mu, kappa, T, M, realization);
    fprintf('         Generate VM noise (bar_gamma_m): 60%%\n');

    bar_gamma_m = function_vmrand(mu, kappa, T, M, realization);
    fprintf('         Generate PE: 80%%\n');
    fprintf('         Generate PE: 100%%\n');

    %% SE Calculation
    disp('------------Step: SE calculation.------------');
    b = 0.2;
    c = 0.43*pi;
    phi_L = -pi/2+c;
    phi_U = pi/2+c;

    %% SE_L MC
    disp('------------Step: SE_L MC.------------');

    beta_SE_L = function_beta(phi_L+bar_delta_m+bar_gamma_m, a, b, c);

    beta_exp_phi_m_SE_L = beta_SE_L .* exp_phi_m .* exp(-1i*gamma_m) .* exp(-1i*delta_m);

    signal_all_realization = zeros(realization, T);

    for aa = 1 : realization
        fprintf('\rSE_L MC progress: %3.0f%%', (aa-1)/realization*100);
        
        [signal_direct] = function_h_d(A_0, h_d_NLoS(:, aa), tau_0, K_1, K_2);
        
        [signal_cascaded] = function_h_RIS(d_AP_RIS, d_RIS_User, T, M, g_m_LoS, h_m_LoS, ...
            beta_exp_phi_m_SE_L(:, :, aa), g_m_NLoS(:, :, aa), h_m_NLoS(:, :, aa),...
            K_3, K_4, G_1, G_2, projected_aperture_RIS_to_AP, projected_aperture_RIS_to_User);
        
        signal_all_realization(aa, :) = signal_cascaded + signal_direct;
        
    end
    fprintf('\rSE_L MC progress: 100%%\n');

    [SE_L_MC, ~] = function_compute_SE(signal_all_realization, Pt_W, sigma_square_W, kappa_t, kappa_r);
    SE_L_MC_values(a_idx) = mean(SE_L_MC); 

    %% SE_U MC
    disp('------------Step: SE_U MC.------------');

    beta_SE_U = function_beta(phi_U+bar_delta_m+bar_gamma_m, a, b, c);

    beta_exp_phi_m_SE_U = beta_SE_U .* exp_phi_m .* exp(-1i*gamma_m) .* exp(-1i*delta_m);

    signal_all_realization = zeros(realization, T);

    for aa = 1 : realization
        fprintf('\rSE_U MC progress: %3.0f%%', (aa-1)/realization*100);
        
        [signal_direct] = function_h_d(A_0, h_d_NLoS(:, aa), tau_0, K_1, K_2);
        
        [signal_cascaded] = function_h_RIS(d_AP_RIS, d_RIS_User, T, M, g_m_LoS, h_m_LoS, ...
            beta_exp_phi_m_SE_U(:, :, aa), g_m_NLoS(:, :, aa), h_m_NLoS(:, :, aa),...
            K_3, K_4, G_1, G_2, projected_aperture_RIS_to_AP, projected_aperture_RIS_to_User);
        
        signal_all_realization(aa, :) = signal_cascaded + signal_direct;
        
    end
    fprintf('\rSE_U MC progress: 100%%\n');

    [SE_U_MC, ~] = function_compute_SE(signal_all_realization, Pt_W, sigma_square_W, kappa_t, kappa_r);
    SE_U_MC_values(a_idx) = mean(SE_U_MC); 

    %% SE MC
    disp('------------Step: SE MC.------------');
    beta_SE = function_beta(phi_m+bar_delta_m+bar_gamma_m, a, b, c);
    beta_exp_phi_m_SE = beta_SE .* exp_phi_m .* exp(-1i*gamma_m) .* exp(-1i*delta_m);
    signal_all_realization = zeros(realization, T);

    for aa = 1 : realization

        percent = (aa-1)/realization*100;
        if percent >= progress_steps(next_step_index)
            fprintf('\rSE MC progress: %3.0f%%', progress_steps(next_step_index));
            if next_step_index < length(progress_steps)
                next_step_index = next_step_index + 1;
            end
        end

        [signal_direct] = function_h_d(A_0, h_d_NLoS(:, aa), tau_0, K_1, K_2);

        [signal_cascaded] = function_h_RIS(d_AP_RIS, d_RIS_User, T, M, g_m_LoS, h_m_LoS, ...
            beta_exp_phi_m_SE(:, :, aa), g_m_NLoS(:, :, aa), h_m_NLoS(:, :, aa),...
            K_3, K_4, G_1, G_2, projected_aperture_RIS_to_AP, projected_aperture_RIS_to_User);

        signal_all_realization(aa, :) = signal_cascaded + signal_direct;
    end

    fprintf('\rSE MC progress: 100%%\n');
    [SE_MC, ~] = function_compute_SE(signal_all_realization, Pt_W, sigma_square_W, kappa_t, kappa_r);
    SE_MC_values(a_idx) = mean(SE_MC); 

end

%% Plotting
close all;
warning off;

disp('------------Step: Plotting the simulation results.------------');

figure(1);
hold on; box on; grid on;

plot(a_values, SE_U_MC_values, ':', 'LineWidth', 3.5, 'Color', [244, 177, 131] ./ 255, 'DisplayName', '$\mathsf{SE}_\mathsf{U}$');
plot(a_values, SE_MC_values, '-.', 'LineWidth', 3.5, 'Color', [128, 116, 200] ./ 255, 'DisplayName', '$\mathsf{SE}$');
plot(a_values, SE_L_MC_values, '-', 'LineWidth', 3.5, 'Color', [169, 209, 142] ./ 255, 'DisplayName', '$\mathsf{SE}_\mathsf{L}$');

xlabel('$a\in[1, 5]$','Interpreter','LaTex');
ylabel('$\mathsf{SE}$ (bit/s/Hz)','Interpreter','LaTex');
legend('Interpreter','LaTex','Location','SouthWest');
xticks(a_values); 

ax = gca;
ax.TickLabelInterpreter = 'latex';   
set(gca,'FontSize',16, 'LineWidth', 1.2);

%% Timer ends
toc;