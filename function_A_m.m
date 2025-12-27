function [A_m_sum, A_m_matrix]=function_A_m(d_AP_RIS, d_RIS_User, T, M, projected_aperture_RIS_to_AP, projected_aperture_RIS_to_User)
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

global lambda_c;


G_RIS_to_AP = ((4*pi) / lambda_c^2) * projected_aperture_RIS_to_AP; % size(projected_aperture_RIS_to_AP) = M x T, projected_aperture_RIS_to_AP = cos(-)*d_x*d_y

G_RIS_to_User = ((4*pi) / lambda_c^2) * projected_aperture_RIS_to_User; % 同上


G_AP_to_RIS = ones(M, T);

A_AP_to_RIS = lambda_c/(4*pi) .* sqrt(G_AP_to_RIS .* G_RIS_to_AP) ./ d_AP_RIS';
% size(A_AP_to_RIS) = 1 x M_total_elements_number

G_User_to_RIS = ones(M, T); % The User is isotropic.

A_RIS_to_User = lambda_c/(4*pi) * sqrt(G_User_to_RIS.*G_RIS_to_User) ./ d_RIS_User';
% size(G_User_to_RIS) = T_total_time_slot x M_total_elements_number
% size(G_User_to_RIS) = T_total_time_slot x M_total_elements_number
% size(d_RIS_User) = T_total_time_slot x M_total_elements_number.
% size(A_RIS_User) = T_total_time_slot x M_total_elements_number;

A_m_sum = sum(A_AP_to_RIS .* A_RIS_to_User, 1)';
% size(A_AP_to_RIS) = 1 x M_total_elements_number;
% size(A_RIS_to_User) = T_total_time_slot x M_total_elements_number;
% size(A_AP_to_RIS .* A_RIS_to_User) = T_total_time_slot x M_total_elements_number;
% size(A_mn_sum) = T_total_time_slot x 1;

A_m_matrix = (A_AP_to_RIS .* A_RIS_to_User)';
% size(A_AP_to_RIS) = 1 x M_total_elements_number
% size(A_RIS_to_User) = T_total_time_slot x M_total_elements_number;
% size(A_mn_matrix) = size(A_AP_to_RIS .* A_RIS_to_User) = T_total_time_slot x M_total_elements_number;

end
