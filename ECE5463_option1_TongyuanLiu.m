%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECE 5463 --> Final Project --> Option 1 
% Tongyuan Liu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Basic parameters
L1 = 10; L2 = 10; L3 = 10;
m1 = 1; m2 = 1; m3 = 1;
g = 9.8;

% PD controler
zeta = 1.0; % from Lecture 17 Page 4 -- Critically damped
wn1 = 20; wn2 = 10; wn3 = 5; % can be changed
Kp = diag([wn1^2, wn2^2 wn3^2]);
Kv = diag([2*zeta*wn1, 2*zeta*wn2 2*zeta*wn3]);
               
% % Position parameter input
% % source from Youtube "ENGR 1315 - MATLAB - inputdlg to Create Kinematics Calculator"
% % https://www.youtube.com/watch?v=61HeXe64A6Y
% prompt = { ...
%     'Enter PICK 1 as "x y":', ...
%     'Enter PLACE 1 as "x y":', ...
%     'Enter PICK 2 as "x y":', ...
%     'Enter PLACE 2 as "x y":' ...
% };
% dlgtitle = 'Input Parameters';
% dims = [1 35]; % input field width (characters)
% definput = {'', '', '', '', ''}; % default values
% input_value = inputdlg(prompt, dlgtitle, dims, definput);
% 
% if isempty(input_value)
% error('User cancelled.');
% end
% 
% pick1_xy = str2num(input_value{1}); % user enters "10 20"
% P_pick_1_xy = pick1_xy(:); % force column vector [x; y]
% place1_xy = str2num(input_value{2});
% P_place_1_xy = place1_xy(:);
% pick2_xy = str2num(input_value{3}); % user enters "10 20"
% P_pick_2_xy = pick2_xy(:); % force column vector [x; y]
% place2_xy = str2num(input_value{4});
% P_place_2_xy = place2_xy(:);

P_pick_1_xy = [15; 0]; P_place_1_xy = [20; 0];
P_pick_2_xy = [25; 0]; P_place_2_xy = [30; 0];
[Pick_1_th1, Pick_1_th2, Pick_1_th3] = IK_component(P_pick_1_xy, L1, L2, L3);
[Place_1_th1, Place_1_th2, Place_1_th3] = IK_component(P_place_1_xy, L1, L2, L3);
[Pick_2_th1, Pick_2_th2, Pick_2_th3] = IK_component(P_pick_2_xy, L1, L2, L3);
[Place_2_th1, Place_2_th2, Place_2_th3] = IK_component(P_place_2_xy, L1, L2, L3);

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

params_eom = struct('L1',L1,'L2',L2,'L3',L3,...
                    'm1',m1,'m2',m2,'m3',m3,'g',g, ...
                    'Kp',Kp,'Kv',Kv, ...
                    'qd_fun',qd_fun,'dqd_fun',dqd_fun,'ddqd_fun',ddqd_fun);

% initial position,code style from my previous code HW4_Q3 and PA3
theta1 = 0; theta2 = 0; theta3 = 0;
dtheta1 = 0; dtheta2 = 0; dtheta3 = 0;
ang = [theta1; theta2; theta3; dtheta1; dtheta2; dtheta3];

% Plot the solution for t=0 to T8
t = [0, T8];
ode45_function_component = @(t,x)closedloop_3R(t,x,params_eom);
[t, x] = ode45(ode45_function_component, t, ang);

% Animated Ploting
% Setting Workspace and plot
figure;
hold on; grid on;
centers = [0, 0];
xlim([0,30]);
ylim([0,30]);
xlabel('X(m)'); ylabel('Y(m)');
title('Final Project Pick-and-Place 3R Planar Arm')


% Inverse Kinematic solution(written by hand, source from HW4 solution)
function [theta1, theta2, theta3] = IK_component(target_xy, L1, L2, L3)
    phi = -deg2rad(45); % constant value, can be changed in code
    wx = target_xy(1) - L3*cos(phi);
    wy = target_xy(2) - L3*sin(phi);
 
    % source from Lecture 11 Page 2
    beta_temp = (L1^2 + L2^2 - wx^2 - wy^2)/(2*L1*L2); % if not doing this, the simulation would be stuck and stopped
    beta_temp = max(min(beta_temp, 1), -1); % stuck between[-1, 1]
    beta = acos(beta_temp);
    alpha_temp = (wx^2 + wy^2 + L1^2 - L2^2)/(2*L1*sqrt(wx^2 + wy^2)); 
    alpha_temp = max(min(alpha_temp, 1), -1);
    alpha = acos(alpha_temp);
    gamma = atan2(wy, wx);
    theta1 = gamma - alpha;
    theta2 = pi - beta;
    theta3 = phi - theta1 - theta2;
    
%     r2 = wx^2 + wy^2;
%     c2 = max(min((r2 - L1^2 - L2^2)/(2*L1*L2), 1), -1);  % cos(theta2)
%     s2 = -sqrt(max(0, 1 - c2^2));  % sin(theta2)
%     theta2 = atan2(s2, c2);
% 
%     k1 = L1 + L2*c2;
%     k2 = L2*s2;
%     theta1 = atan2(wy, wx) - atan2(k2, k1);
%     theta3 = phi - theta1 - theta2;
%     
%     theta1 = mod(theta1 + pi, 2*pi) - pi;
%     theta2 = mod(theta2 + pi, 2*pi) - pi;
%     theta3 = mod(theta3 + pi, 2*pi) - pi;
    
end

% Use element from EoM to build PD controler, equation source from Lecture17 Page 4 and Page 9
% Also use ode for later plot
function dx = closedloop_3R(t,x,params)
    % origial value
    q = x(1:3); dq = x(4:6); 
    
    % expected value
    qd = params.qd_fun(t);
    dqd = params.dqd_fun(t);
    ddqd = params.ddqd_fun(t);
    if numel(qd)~=3 || numel(dqd)~=3 || numel(ddqd)~=3
        error('qd_fun/dqd_fun/ddqd_fun must each return a 3x1 vector.');
    end
    if any(~isfinite([qd; dqd; ddqd]))
        error('qd/dqd/ddqd contain NaN/Inf. Check your trajectory functions and inputs.');
    end
    
    % error calculation
    err = qd - q; derr = dqd - dq;
    
    % calculate EoM element
    [M,C,G] = EoM_3R_terms([q,dq],...
                            params.L1,params.L2,params.L3,...
                            params.m1,params.m2,params.m3,...
                            params.g);
    if ~isequal(size(params.Kp), [3 3]) || ~isequal(size(params.Kv), [3 3])
        error('Kp and Kv must be 3x3 matrices.');
    end
    
    % calculate the velocity, use velocity to compute tau
    % source from Lecture17 Page 9
    V = ddqd + params.Kv*derr + params.Kp*err;
    tau = M*V + C + G;
    
    ddq = M\(tau - C - G);
    dx = [dq; ddq];
end

% Use EoM to calculate theta derivative, method from Lecture 11
% source from Youtube "Equations of Motion for the Multi Degree of Freedom (MDOF) Problem Using LaGrange's Equations"
% https://www.youtube.com/watch?v=uAKD5CGZuSs
function [M,C,G] = EoM_3R_terms(x,L1,L2,L3,m1,m2,m3,g)
    % x = [q1; q2; q3; dq1; dq2; dq3]
    q1 = x(1); q2 = x(2); q3 = x(3);
    dq1 = x(4); dq2 = x(5); dq3 = x(6);
    
    % calculate Mass matrix
    M = zeros(3,3);
    a1 = (m1*L1^2)/3 + m2*L1^2 + m3*L1^2;  
    a2 = (m2*L2^2)/3 + m3*L2^2;   
    a3 = (m3*L3^2)/3; 
    b12 = m2*L1*(L2/2) + m3*L1*L2;              
    b13 = m3*L1*(L3/2);                           
    b23 = m3*L2*(L3/2);                       
    M(1,1) = a1 + a2 + a3 + 2*b12*cos(q2) + 2*b13*cos(q2+q3) + 2*b23*cos(q3);
    M(1,2) = a2 + a3 + b12*cos(q2) + b13*cos(q2+q3) + b23*cos(q3);
    M(1,3) = a3 + b13*cos(q2+q3) + b23*cos(q3);
    M(2,1) = M(1,2);
    M(2,2) = a2 + a3 + 2*b23*cos(q3);
    M(2,3) = a3 + b23*cos(q3);
    M(3,1) = M(1,3);
    M(3,2) = M(2,3);
    M(3,3) = a3;
    
    % calculate velocity product terms
    C = zeros(3,1);
    h1 = -b12*sin(q2);            
    h2 = -b13*sin(q2+q3);        
    h3 = -b23*sin(q3);    
    C(1,1) = (2*h1*dq1*dq2 + h1*dq2^2)...
            + (2*h2*dq1*(dq2+dq3) + h2*(dq2+dq3)^2) ...
            + (2*h3*dq1*dq3 + h3*dq3^2);
    C(2,1) = (-h1*dq1^2) ...
            + (2*(-b23*sin(q3)) * dq2*dq3 + (-b23*sin(q3))*dq3^2) ...
            + (-b13*sin(q2+q3))*(dq1^2 + 2*dq1*dq2 + (dq2)^2);
    C(3,1) = (-h2)*(dq1^2 + 2*dq1*dq2 + (dq2)^2) ...
            + (-h3)*(dq1^2 + 2*dq1*dq3 + 2*dq2*dq3 + dq3^2);
    
    % calculate gravity terms
    G = zeros(3,1);
    G(1,1) = g*(m1*(L1/2)*cos(q1) ...
           + m2*(L1*cos(q1) + (L2/2)*cos(q1+q2)) ...
           + m3*(L1*cos(q1) + L2*cos(q1+q2) + (L3/2)*cos(q1+q2+q3)));
    G(2,1) = g*(m2*((L2/2)*cos(q1+q2)) ...
           + m3*(L2*cos(q1+q2) + (L3/2)*cos(q1+q2+q3)));
    G(3,1) = g*(m3*((L3/2)*cos(q1+q2+q3)));
end

% Calculate theta, dtheta and ddtheta
function [q, dq, ddq] = P2P_trajectory(t, t0, tf, q0, qf) 
    % For both initial and end, q,dq and ddq have no different change
    if t <= t0
        q = q0; dq = zeros(3,1); ddq = zeros(3,1);
    elseif t >= tf
        q = qf; dq = zeros(3,1); ddq = zeros(3,1);
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
    output = q(:);
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
        dq = [0;0;0]; % on the point "Pick 1", and stay, no movement
    elseif t <= T3
        [~,dq,~] = P2P_trajectory(t, T2, T3, P_pick_1, P_place_1);
    elseif t <= T4
        dq = [0;0;0]; % on the point "Place 1", and stay, no movement
    elseif t <= T5
        [~,dq,~] = P2P_trajectory(t, T4, T5, P_place_1, P_pick_2);
    elseif t <= T6
        dq = [0;0;0]; % on the point "Pick 2", and stay, no movement
    elseif t <= T7
        [~,dq,~] = P2P_trajectory(t, T6, T7, P_pick_2, P_place_2);
    else
        dq = [0;0;0]; % on the point "Place 2", and stay, no movement
    end
    output = dq(:);
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
        ddq = [0;0;0]; % on the point "Pick 1", and stay, no movement
    elseif t <= T3
        [~,~,ddq] = P2P_trajectory(t, T2, T3, P_pick_1, P_place_1);
    elseif t <= T4
        ddq = [0;0;0]; % on the point "Place 1", and stay, no movement
    elseif t <= T5
        [~,~,ddq] = P2P_trajectory(t, T4, T5, P_place_1, P_pick_2);
    elseif t <= T6
        ddq = [0;0;0]; % on the point "Pick 2", and stay, no movement
    elseif t <= T7
        [~,~,ddq] = P2P_trajectory(t, T6, T7, P_pick_2, P_place_2);
    else
        ddq = [0;0;0]; % on the point "Place 2", and stay, no movement
    end
    output = ddq(:);
end



