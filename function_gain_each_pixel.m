function channel_gain_each_pixel = function_gain_each_pixel(D_AP_or_User, D_m, dx, dy)
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

t_1 = D_m(1) - dx/2 - D_AP_or_User(1);
t_2 = D_m(1) + dx/2 - D_AP_or_User(1);
t_3 = D_m(2) - dy/2 - D_AP_or_User(2);
t_4 = D_m(2) +dy/2 - D_AP_or_User(2);
z = D_AP_or_User(3);

Q = @(s_1, s_2, z) 1/z * atan(s_1*s_2/(z*sqrt(s_1^2+s_2^2+z^2)));

channel_gain_each_pixel = Q(t_2, t_4, z) - Q(t_1, t_4, z) - Q(t_2, t_3, z) + Q(t_1, t_3, z);

end 
