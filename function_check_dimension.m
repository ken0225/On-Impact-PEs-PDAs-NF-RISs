function [coordinate_output] = function_check_dimension(coordinate_input)
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

% Obtain the dimension (size) information of the input coordinate matrix.
size_coordinate_input = size(coordinate_input);

% Check if the input is a matrix
if ismatrix(coordinate_input)
    
    % Check if the size of the input is (X x 3).
    if size_coordinate_input(2) == 3
        
        % If the size of the input is (X x 3), output is the same.
        coordinate_output = coordinate_input;
        
        % If not, check if the size of the input is (3 x X).
    elseif size_coordinate_input(1) == 3
        
        % If the size of the input is (3 x X), transform it to (input.').
        coordinate_output = coordinate_input.';
        
    else
        
        % If the size of the input is (X x X), error appears.
        error('The input is NOT an Euclidean vector.')
        
    end
    
else
    
    % If the input is not a matrix, error appears.
    error('Only matrix supported.')
    
end

end
