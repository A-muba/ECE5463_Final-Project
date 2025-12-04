%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECE 5463 --> Final Project --> Option 1 
% Tongyuan Liu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Basic parameters
L1 = 1; L2 = 1; L3 = 1;
m1 = 1; m2 = 1; m3 = 1;
g = 9.8;
zeta = 1.0; % from Lecture 17 Page 4 -- Critically damped
               
% Position parameter input
% source from Youtube "ENGR 1315 - MATLAB - inputdlg to Create Kinematics Calculator"
% https://www.youtube.com/watch?v=61HeXe64A6Y
prompt = { ...
    'Enter phi (deg) (phi is the sum of 3 theta):', ...
    'Enter PICK 1 as "x y":', ...
    'Enter PLACE 1 as "x y":', ...
    'Enter PICK 2 as "x y":', ...
    'Enter PLACE 2 as "x y":' ...
};
dlgtitle = 'Input Parameters';
dims = [1 35]; % input field width (characters)
definput = {'', '', '', '', ''}; % default values
input_value = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(input_value)
error('User cancelled.');
end

phi_deg = str2double(input_value{1});
pick1_xy = str2double(input_value{2}); % user enters "10 20"
P_pick_1 = pick1_xy(:); % force column vector [x; y]
place1_xy = str2double(input_value{3});
P_place_1 = place1_xy(:);
pick2_xy = str2double(input_value{4}); % user enters "10 20"
P_pick_2 = pick2_xy(:); % force column vector [x; y]
place2_xy = str2double(input_value{5});
P_place_2 = place2_xy(:);

% Inverse Kinematic solution(written by hand, source from HW4 solution)
function [theta1, theta2, theta3] = IK_component(phi_deg, target_xy)
    Px = target_xy(1); Py = target_xy(2);
    % source from Lecture 11 Page 2
    beta_temp = (L1^2 + L2^2 - Px^2 - Py^2)/(2*L1*L2); % if not doing this, the simulation would be stuck and stopped
    beta_temp = max(min(beta_temp, 1), -1); % stuck between[-1, 1]
    beta = acos(beta_temp);
    alpha_temp = (Px^2 + Py^2 + L1^2 - L2^2)/(2*L1*sqrt(Px^2 + Py^2)); 
    alpha_temp = max(min(alpha_temp, 1), -1);
    alpha = acos(alpha_temp);
    gamma = atan2(Py, Px);
    theta1 = gamma - alpha;
    theta2 = pi - beta;
    theta3 = phi_deg - theta1 - theta2;
end









