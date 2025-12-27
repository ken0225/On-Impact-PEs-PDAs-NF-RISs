%% Program Information
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

function vn = function_refcoefficient(omega,Cn,Rn)
% Compute the reflection coefficient using the circuit from (3) in
% "Intelligent Reflecting Surface: Practical Phase Shift Model and
% Beamforming Optimization" by Samith Abeywickrama, Rui Zhang, Chau Yuen.

% Input: 
% omega, i.e., carrier frequency
% Cn, i.e., effective capacitance
% Rn, i.e., effective resistance, Rn = [1 1.5 2 2.5]

% Output: vn, i.e., reflection coefficient

L1 = 2.5e-9; %Bottom layer inductance
L2 = 0.7e-9; %Top layer inductance
Z0 = 377; %Impedance of free space

%Compute the impedance of the surface
Zn = (1i*omega*L1*(1i*omega*L2+1./(1i*omega*Cn)+Rn))./(1i*omega*L1+(1i*omega*L2+1./(1i*omega*Cn)+Rn));

%Compute reflection coefficient
vn = (Zn-Z0)./(Zn+Z0);

end