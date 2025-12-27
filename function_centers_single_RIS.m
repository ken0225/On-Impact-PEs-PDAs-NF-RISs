function [output_centers_RIS]=function_centers_single_RIS(X, Y, d_x, d_y, h_RIS)
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

X_interval = function_interval(X);

Y_interval = function_interval(Y);

X_range = (X_interval(1) : X_interval(2));

Y_range = (Y_interval(1) : Y_interval(2));

for aa = 1 : length(X_range)
    
    g_X(aa) = X_range(aa)*d_x - 0.5*d_x*mod(X+1, 2);
    
end

for bb = 1 : length(Y_range)
    
    g_Y(bb) = Y_range(bb)*d_y - 0.5*d_y*mod(Y+1, 2);
    
end

temp1_output_centers_RIS = function_cartprod(g_X, g_Y);

temp1_output_centers_RIS(:, end+1) = 0;

temp2_output_centers_RIS = sortrows(temp1_output_centers_RIS, 1);

temp2_output_centers_RIS(:, 2) = temp2_output_centers_RIS(:, 2) + h_RIS;

[output_centers_RIS] = function_check_dimension(temp2_output_centers_RIS);

end

function [output_interval]=function_interval(S)
% Author: Ke(Ken)WANG from Macao Polytechnic University.
% Email: ke.wang@mpu.edu.mo, kewang0225@gmail.com.
% Update infomation: v0.1(2020/Mar), v0.2(2021/Sep), v0.3(2023/Mar/26).
%
% Suppose we have an RIS with M x N elements, this function aims to
% calculate the smallest and biggest indices of M and N.
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% original article.
%
% ===========================================
%
% Example:
%
% [output_interval]=function_interval(8)
%
% output_interval =
%
%     -3     4
%
% [output_interval]=function_interval(11.2)
%
% the edges of interval are not integer!
%
% [output_interval]=function_interval(11)
%
% output_interval =
%
%     -5     5
%
% ===========================================

S_min = mod(S+1, 2) - floor(S/2);

S_max = floor(S/2);

assert(S_min == fix(S_min) & S_max == fix(S_max), 'the edges of interval are not integer!');

output_interval = [fix(S_min), fix(S_max)];

end
