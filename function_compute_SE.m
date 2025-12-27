function [mean_SE_MC, mean_SNR_MC]=function_compute_SE(all_realization, Pt_transmit_power, sigma_square_noise_power_W, kappa_t, kappa_r)
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

received_power_all_realization = Pt_transmit_power .* mean((abs(all_realization)) .^ 2);

mean_SNR_MC = received_power_all_realization ./ ...
    (kappa_t .* received_power_all_realization + kappa_r .* received_power_all_realization + sigma_square_noise_power_W);

mean_SE_MC = log2(1 + mean_SNR_MC);

end