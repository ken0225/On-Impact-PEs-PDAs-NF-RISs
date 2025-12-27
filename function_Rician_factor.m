function [kappa_d, rho_m, kappa_m]=function_Rician_factor(d_0, d_AP_RIS, d_RIS_User, varrho, iota)
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

kappa_d = db2pow(varrho - iota.*d_0);

rho_m = db2pow(varrho - iota.*d_AP_RIS);

kappa_m = db2pow(varrho - iota.*d_RIS_User);

end
