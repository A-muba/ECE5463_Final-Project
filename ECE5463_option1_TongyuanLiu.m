%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECE 5463 --> Final Project --> Option 1 
% Tongyuan Liu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Basic parameters
L1 = 1; L2 = 1; L3 = 1;
m1 = 1; m2 = 1; m3 = 1;
g = 9.8;

% PD controler
zeta = 1.0; % from Lecture 17 Page 4 -- Critically damped
wn1 = 20; wn2 = 10; % can be changed
Kp = diag([wn1^2, wn2^2]);
Kv = diag([2*zeta*wn1, 2*zeta*wn2]);
               
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
P_pick_1_xy = pick1_xy(:); % force column vector [x; y]
place1_xy = str2double(input_value{3});
P_place_1_xy = place1_xy(:);
pick2_xy = str2double(input_value{4}); % user enters "10 20"
P_pick_2_xy = pick2_xy(:); % force column vector [x; y]
place2_xy = str2double(input_value{5});
P_place_2_xy = place2_xy(:);

[Pick_1_th1, Pick_1_th2, Pick_1_th3] = IK_component(phi_deg, P_pick_1_xy);
[Place_1_th1, Place_1_th2, Place_1_th3] = IK_component(phi_deg, P_place_1_xy);
[Pick_2_th1, Pick_2_th2, Pick_2_th3] = IK_component(phi_deg, P_pick_2_xy);
[Place_2_th1, Place_2_th2, Place_2_th3] = IK_component(phi_deg, P_place_2_xy);

% position initial
P_org_rad = [0; 0; 0];                     
P_pick_1_rad = [deg2rad(Pick_1_th1); deg2rad(Pick_1_th2); deg2rad(Pick_1_th3)];
P_place_1_rad = [deg2rad(Place_1_th1); deg2rad(Place_1_th2); deg2rad(Place_1_th3)];
P_pick_2_rad = [deg2rad(Pick_2_th1); deg2rad(Pick_2_th2); deg2rad(Pick_2_th3)];
P_place_2_rad = [deg2rad(Place_2_th1); deg2rad(Place_2_th2); deg2rad(Place_2_th3)];

% Position holding time
% P_org->A, Pick1->B, Place1->C, Pick2->D, Place2->E
T_A2B = 1.5;
T_holdB = 1;
T_B2C = 1.5; 
T_holdC = 1;
T_C2D = 1.5;
T_holdD = 1;
T_D2E = 1.5;
T_holdE = 1;

% total lasting time, for polt
T1 = T_A2B; 
T2 = T_A2B + T_holdB;
T3 = T_A2B + T_holdB + T_B2C;
T4 = T_A2B + T_holdB + T_B2C + T_holdC;
T5 = T_A2B + T_holdB + T_B2C + T_holdC + T_C2D;
T6 = T_A2B + T_holdB + T_B2C + T_holdC + T_C2D + T_holdD;
T7 = T_A2B + T_holdB + T_B2C + T_holdC + T_C2D + T_holdD + T_D2E;
T8 = T_A2B + T_holdB + T_B2C + T_holdC + T_C2D + T_holdD + T_D2E + T_holdE;

% Get the theta parameter for later calculation
% the frameware of the algorithm is from my own work PA3 with different parameter 
params_destination = struct('P_org',P_org_rad,...
                            'P_pick_1',P_pick_1_rad,'P_place_1',P_place_1_rad,...
                            'P_pick_2',P_pick_2_rad,'P_place_2',P_place_2_rad,...
                            'T1',T1,'T2',T2,'T3',T3,'T4',T4,'T5',T5,'T6',T6,...
                            'T7',T7);
qd_fun = @(t) qd_piece(t, params_destination);
dqd_fun = @(t) dqd_piece(t, params_destination);
ddqd_fun = @(t) ddqd_piece(t, params_destination);

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

% Calculate theta, dtheta and ddtheta
function [q, dq, ddq] = P2P_trajectory(t, t0, tf, q0, qf) 
    % For both initial and end, q,dq and ddq have no different change
    if t <= t0
        q = q0; dq = zeros(2,1); ddq = zeros(2,1);
    elseif t >= tf
        q = qf; dq = zeros(2,1); ddq = zeros(2,1);
    % for the period, I used the method in the text book, also shown in the
    % video. Source: Modern Robotics, Chapters 9.1 and 9.2: Point-to-Point Trajectories (Part 2 of 2)
    % https://www.youtube.com/watch?v=0ZqeBEa_MWo
    else
        tau = (t - t0)/(tf - t0);
        s = 10*tau^3 - 15*tau^4 + 6*tau^5;
        ds= (30*tau^2 - 60*tau^3 + 30*tau^4)/(tf - t0); % derivative of s
        dds = (60*tau - 180*tau^2 + 120*tau^3)/((tf - t0)^2); % derivative of ds
        q = q0 + (qf - q0)*s;
        dq = (qf - q0)*ds; % derivative of q
        ddq = (qf - q0)*dds; % derivative of dq
    end
end

% for q, dq and ddq, set 3 function to illustate
% the frameware of the algorithm is from my own work PA3 with different parameter 
function output = qd_piece(t, params_destination)
    P_org = params_destination.P_org;
    P_pick_1 = params_destination.P_pick_1;
    P_place_1 = params_destination.P_place_1;
    P_pick_2 = params_destination.P_pick_2;
    P_place_2 = params_destination.P_place_2;
    T1 = params_destination.T1; T2 = params_destination.T2;
    T3 = params_destination.T3; T4 = params_destination.T4;
    T5 = params_destination.T5; T6 = params_destination.T6;
    T7 = params_destination.T7;
    
    if t <= T1
        [q,~,~] = P2P_trajectory(t, 0, T1, P_org, P_pick_1);
    elseif t <= T2
        q = P_pick_1; % on the point "Pick 1", and stay
    elseif t <= T3
        [q,~,~] = P2P_trajectory(t, T2, T3, P_pick_1, P_place_1);
    elseif t <= T4
        q = P_place_1; % on the point "Place 1", and stay
    elseif t <= T5
        [q,~,~] = P2P_trajectory(t, T4, T5, P_place_1, P_pick_2);
    elseif t <= T6
        q = P_pick_2; % on the point "Pick 2", and stay
    elseif t <= T7
        [q,~,~] = P2P_trajectory(t, T6, T7, P_pick_2, P_place_2);
    else
        q = P_place_2; % on the point "Place 2", and stay
    end
    output = q;
end

% the frameware of the algorithm is from my own work PA3 with different parameter 
function output = dqd_piece(t, params_destination)
    P_org = params_destination.P_org;
    P_pick_1 = params_destination.P_pick_1;
    P_place_1 = params_destination.P_place_1;
    P_pick_2 = params_destination.P_pick_2;
    P_place_2 = params_destination.P_place_2;
    T1 = params_destination.T1; T2 = params_destination.T2;
    T3 = params_destination.T3; T4 = params_destination.T4;
    T5 = params_destination.T5; T6 = params_destination.T6;
    T7 = params_destination.T7;
    
    if t <= T1
        [~,dq,~] = P2P_trajectory(t, 0, T1, P_org, P_pick_1);
    elseif t <= T2
        dq = [0;0]; % on the point "Pick 1", and stay, no movement
    elseif t <= T3
        [~,dq,~] = P2P_trajectory(t, T2, T3, P_pick_1, P_place_1);
    elseif t <= T4
        dq = [0;0]; % on the point "Place 1", and stay, no movement
    elseif t <= T5
        [~,dq,~] = P2P_trajectory(t, T4, T5, P_place_1, P_pick_2);
    elseif t <= T6
        dq = [0;0]; % on the point "Pick 2", and stay, no movement
    elseif t <= T7
        [~,dq,~] = P2P_trajectory(t, T6, T7, P_pick_2, P_place_2);
    else
        dq = [0;0]; % on the point "Place 2", and stay, no movement
    end
    output = dq;
end

% the frameware of the algorithm is from my own work PA3 with different parameter 
function output = ddqd_piece(t, params_destination)
    P_org = params_destination.P_org;
    P_pick_1 = params_destination.P_pick_1;
    P_place_1 = params_destination.P_place_1;
    P_pick_2 = params_destination.P_pick_2;
    P_place_2 = params_destination.P_place_2;
    T1 = params_destination.T1; T2 = params_destination.T2;
    T3 = params_destination.T3; T4 = params_destination.T4;
    T5 = params_destination.T5; T6 = params_destination.T6;
    T7 = params_destination.T7;
    
    if t <= T1
        [~,~,ddq] = P2P_trajectory(t, 0, T1, P_org, P_pick_1);
    elseif t <= T2
        ddq = [0;0]; % on the point "Pick 1", and stay, no movement
    elseif t <= T3
        [~,~,ddq] = P2P_trajectory(t, T2, T3, P_pick_1, P_place_1);
    elseif t <= T4
        ddq = [0;0]; % on the point "Place 1", and stay, no movement
    elseif t <= T5
        [~,~,ddq] = P2P_trajectory(t, T4, T5, P_place_1, P_pick_2);
    elseif t <= T6
        ddq = [0;0]; % on the point "Pick 2", and stay, no movement
    elseif t <= T7
        [~,~,ddq] = P2P_trajectory(t, T6, T7, P_pick_2, P_place_2);
    else
        ddq = [0;0]; % on the point "Place 2", and stay, no movement
    end
    output = ddq;
end



