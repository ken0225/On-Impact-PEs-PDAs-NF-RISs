function [output_projected_angle]=function_projected_angle(p1, p2) 
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

    %Note that "(p2(:, end) - p1(:, end)" is the z-axis coordinate 
    %output_projected_angle = acos((p2(:, end)-p1(:, end)) ./ vecnorm((p2-p1).').');
    output_projected_angle = acos((p2(:, 3)-p1(:, 3)) ./ vecnorm(p2-p1,2,2) );
    
    % Convert the column vector to row vector
    output_projected_angle = output_projected_angle.'; 
    
end
