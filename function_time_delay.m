function [tau_0, tau_mn, tau_BS_RIS, tau_RIS_MR, d_0, d_BS_RIS, d_RIS_MR]=function_time_delay(tensor_D_centers_RIS, D_BS, D_MR)
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

global c_0;

%% Compute

% Compute the direct link

d_0 = vecnorm((D_MR - D_BS), 2); % size(d_0) = 1 x 1;

tau_0 = d_0 / c_0; % size(tau_0) = 1 x 1;

% Compute the cascaded link

tau_BS_RIS = [];

tau_RIS_MR = [];

d_BS_RIS = [];

d_RIS_MR = [];

for aa = 1 : size(tensor_D_centers_RIS, 3) % size(trajectory_tensor_centers_RIS, 3) = T_sampling_point.
    
    temp_centers_RIS = tensor_D_centers_RIS(:,:,aa);
    
    for bb = 1 : size(tensor_D_centers_RIS, 1) % size(trajectory_tensor_centers_RIS, 1) = size(temp_centers_RIS, 1) = M_total_elements_number.
        
        % Calculate the propagation time between the mn-th element of the RIS and the BS.
        temp1 = vecnorm((temp_centers_RIS(bb, :) - D_BS));
        
        temp_d_BS_RIS(:, bb) = temp1;
        
        temp_tau_BS_RIS(:, bb) = temp1 / c_0;
        % size(tau_BS_RIS) = 1 x M_total_elements_number.
        
        % Calculate the delay between the MR and the mn-th element of the RIS every time slot.
        temp2 = vecnorm((D_MR - temp_centers_RIS(bb, :)), 2);
        % size(temp2) = T_total_time_slot x 1.
        
        temp_d_RIS_MR(:, bb) = temp2;
        % size(d_RIS_MR) = T_total_time_slot x M_total_elements_number.
        
        temp_tau_RIS_MR(:, bb) = temp2 / c_0;
        % size(tau_RIS_MR) = T_total_time_slot x M_total_elements_number.
        
        temp_tau_mn(:, bb) =  (temp1 + temp2) / c_0;
        % size(tau_mn) = T_total_time_slot x M_total_elements_number.
        
    end
    
    d_BS_RIS(aa,:) = temp_d_BS_RIS;
    
    tau_BS_RIS(aa,:) = temp_tau_BS_RIS;
    
    d_RIS_MR(aa,:) = temp_d_RIS_MR;
    
    tau_RIS_MR(aa,:) = temp_tau_RIS_MR;
    
    tau_mn(aa,:) = temp_tau_mn;
    
end
