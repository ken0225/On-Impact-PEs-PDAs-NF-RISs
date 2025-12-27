function [A_0]=function_A_0(d_0)
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

global c_0 f_c;

lambda_c = c_0/f_c;

%The gain of transmitter, from the direction of the receiver.
G_T_to_R = 1;

%The gain of receiver, from the direction of the transmitter.
G_R_to_T = 1;

A_0 = lambda_c/(4*pi) .* sqrt(G_T_to_R.*G_R_to_T) ./ d_0;

end
