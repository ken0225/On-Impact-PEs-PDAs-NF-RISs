function [h_RIS]=function_h_RIS(d_BS_RIS, d_RIS_MR, T, M, g_mn_LoS, h_mn_LoS, exp_phi_mn,...
    g_mn_NLoS_each, h_mn_NLoS_each, K_3, K_4, G_1, G_2, projected_aperture_RIS_to_BS, projected_aperture_RIS_to_MR)
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

global  lambda_c;

G_RIS_to_BS = ((4*pi) / lambda_c^2) * projected_aperture_RIS_to_BS; % size(projected_aperture_RIS_to_BS) = M x T

G_RIS_to_MR = ((4*pi) / lambda_c^2) * projected_aperture_RIS_to_MR;


%% BS to RIS Link

G_BS_to_RIS = ones(M, T);

A_BS_to_RIS = lambda_c/(4*pi) * sqrt(G_BS_to_RIS .* G_RIS_to_BS) ./ d_BS_RIS'; 
% size(A_BS_to_RIS) = M x T.
% size(d_BS_RIS') = M x T.

g_m = [];

for aa = 1 : T
    
    temp_g_m =  A_BS_to_RIS(:, aa)' .* (G_1(aa) .* g_mn_LoS(aa, :) + G_2(aa) .* g_mn_NLoS_each(:, aa).'); % size(temp_g_m) = 1 x M_total_elements_number.
    % A_BS_to_RIS(:, aa) = M x 1
    % size(g_mn_LoS(1, :)) = 1 x M
    % size(G_1) = 1 x 1.
    % size(g_mn_LoS) = T x M
    % size(exp_BS_RIS) = 1 x M_total_elements_number.
    % size(g_mn_NLoS) = M_total_elements_number x T_total_time_slot x K_realization_number.
    % size(g_mn_NLoS_each) = size(g_mn_NLoS(: ,:, aa)) = M_total_elements_number x T_total_time_slot.
    % size(g_mn_NLoS_each(:, aa).') = 1 x M_total_elements_number.
    
    g_m(aa, :) = temp_g_m;
    % size(temp_g_m) = 1 x M_total_elements_number.
    % size(g_m) = T_total_time_slot x M_total_elements_number.
    
    clear temp_g_m
    
end
%% RIS to MR Link

G_MR_to_RIS = ones(M, T); % The MR is isotropic.


h_m= [];

A_RIS_to_MR = lambda_c/(4*pi) * sqrt(G_MR_to_RIS.*G_RIS_to_MR) ./ d_RIS_MR';
% size(G_MR_to_RIS) = T_total_time_slot x M_total_elements_number.
% size(G_MR_to_RIS) = T_total_time_slot x M_total_elements_number.
% size(A_RIS_MR) = T_total_time_slot x M_total_elements_number.
% size(d_RIS_MR) = T_total_time_slot x M_total_elements_number.

for bb = 1 : T
    
    temp_h_m = A_RIS_to_MR(:, bb)' .* (K_3(bb) .* h_mn_LoS(bb, :) + K_4(bb) .* (h_mn_NLoS_each(:, bb).'));
    % size(A_RIS_MR(bb, :)) = 1 x M_total_elements_number.
    % size(h_mn_NLoS(:, :, bb)) = M_total_elements_number x T_total_time_slot.
    % size(h_mn_NLoS_each(:, bb).') = 1 x M_total_elements_number.
    % size(K_3(bb) .* (exp_RIS_MR(bb, :))) = 1 x M_total_elements_number.
    % size((h_mn_NLoS_each(:, bb).')) = 1 x M_total_elements_number.
    
    h_m(bb, :) = temp_h_m;
    % size(temp_h_m) = 1 x M_total_elements_number.
    % size(h_m) = T_total_time_slot x M_total_elements_number.
    
    clear temp_h_m
    
end

%% Compute h_RIS

h_RIS = [];

for cc = 1 : T
    
    %temp_h_RIS = sum(g_m(cc, :) .* h_m(cc, :) .* exp_phi_mn(cc, :) .* exp(-1i*2*pi*f_c));
    temp_h_RIS = sum(g_m(cc, :) .* h_m(cc, :) .* exp_phi_mn(cc, :));
    % size(g_m) = T_total_time_slot x M_total_elements_number.
    % size(h_m) = T_total_time_slot x M_total_elements_number.
    % size(exp_phi_mn) = T_total_time_slot x M_total_elements_number.
    % size(exp(-1i*2*pi*(f_c+f_mn(cc, :)))) = 1 x M_total_elements_number.
    
    h_RIS(cc, :) = temp_h_RIS;
    % size(temp_h_RIS) = 1 x M_total_elements_number.
    % size(h_RIS) = T_total_time_slot x M_total_elements_number.
    
    clear temp_h_RIS
    
end

end
