clc;clear;

% Import the URDF and define the robotic manipulator
robot = importrobot('arm.urdf');
robot.Gravity = [0; 0; -9.81]; % Gravitational acceleration
robot.DataFormat = 'column';   

% Simulation parameters
dt = 0.01;      % Time step
t_end = 0.5;     % Simulation duration
t = 0:dt:t_end;   % Time vector

% Initial conditions
numJoints = numel(robot.homeConfiguration); % Number of joints
q = zeros(numJoints, 1);   % Joint positions (rad)
qd = [1;2;2;-1];  % Joint velocities, initial values (rad/s)
qdd = zeros(numJoints, 1);

% Define screw axes for joints in the base frame
S_all = [...
    0  0  1  0    0      0;
    0  1  0  -0.1775 0   0;
    0  1  0  -0.4945 0   0;
    0  1  0  -0.6965 0   0]';

% Home Configuartion Matrices
M1_0 = [1 0 0 0;0 1 0 0;0 0 1 0.0545;0 0 0 1];
M2_0 = [1 0 0 0;0 1 0 0;0 0 1 0.1775;0 0 0 1];
M3_0 = [1 0 0 0;0 1 0 0;0 0 1 0.4945;0 0 0 1];
M4_0 = [1 0 0 0;0 1 0 0;0 0 1 0.6965;0 0 0 1];

Ad_M1 = computeAdjoint(inv(M1_0));
Ad_M2 = computeAdjoint(inv(M2_0));
Ad_M3 = computeAdjoint(inv(M3_0));
Ad_M4 = computeAdjoint(inv(M4_0));

A1 = Ad_M1*S_all(:,1);
A2 = Ad_M2*S_all(:,2);
A3 = Ad_M3*S_all(:,3);
A4 = Ad_M4*S_all(:,4);

M0_1 = eye(4);
M1_2 = M1_0\M2_0;
M2_3 = M2_0\M3_0;
M3_4 = M3_0\M4_0;

M2_1 = M2_0\M1_0;
M3_2 = M3_0\M2_0;
M4_3 = M4_0\M3_0;

A_all = blkdiag(A1, A2, A3, A4);

A1_skew = screwToMatrix(A1);
A2_skew = screwToMatrix(A2);
A3_skew = screwToMatrix(A3);
A4_skew = screwToMatrix(A4);

% Define spatial inertias
G_all = blkdiag( ...
    diag([0.00028, 0.00027, 0.000065, 0.19, 0.19, 0.19]), ...
    diag([0.00251, 0.00248, 0.000092, 0.293, 0.293, 0.293]), ...
    diag([0.00079, 0.00077, 0.000069, 0.22, 0.22, 0.22]), ...
    diag([0.00037, 0.00036, 0.0000496, 0.16, 0.16, 0.16])...
);

G1 = G_all(1:6,1:6);
G2 = G_all(7:12,7:12);
G3 = G_all(13:18,13:18);
G4 = G_all(19:24,19:24);

% Base velocity and acceleration
V_base = zeros(6, 1);
dV_base = [0; 0; 0; 0; 0; -9.81];

% Simulation loop
for i = 1:length(t)
    % Step 1: Forward Pass (Compute Velocities and Accelerations)
    V0 = V_base;
    dV0 = dV_base;
   
    T0_1 = M0_1 * expm(A1_skew*q(1));
    T1_2 = M1_2 * expm(A2_skew*q(2));
    T2_3 = M2_3 * expm(A3_skew*q(3));
    T3_4 = M3_4 * expm(A4_skew*q(4));

    T1_0 = expm(-A1_skew*q(1)) * M1_0;
    T2_1 = expm(-A2_skew*q(2)) * M2_1;
    T3_2 = expm(-A3_skew*q(3)) * M3_2;
    T4_3 = expm(-A4_skew*q(4)) * M4_3;
    T5_4 = [1 0 0 0;0 1 0 0;0 0 1 -0.1605;0 0 0 1];
        
    % Spatial velocity and acceleration
    V1 = A1 * qd(1) + computeAdjoint(T0_1) * V0;
    V2 = A2 * qd(2) + computeAdjoint(T1_2) * V1;
    V3 = A3 * qd(3) + computeAdjoint(T2_3) * V2;
    V4 = A4 * qd(4) + computeAdjoint(T3_4) * V3;

    dV1 = A1 * qdd(1) + computeAdjoint(T0_1) * dV0 + ad(V1) * A1 * qd(1);
    dV2 = A2 * qdd(2) + computeAdjoint(T1_2) * dV1 + ad(V2) * A2 * qd(2);
    dV3 = A3 * qdd(3) + computeAdjoint(T2_3) * dV2 + ad(V3) * A3 * qd(3);
    dV4 = A4 * qdd(4) + computeAdjoint(T3_4) * dV3 + ad(V4) * A4 * qd(4);   

    % Step 2: Backward Pass (Compute Wrenches and Torques)
    F5 = zeros(6, 1); % No external wrench at the end-effector
    
    F4 = G4 * dV4 - ad(V4).' * G4 * V4 + computeAdjoint(T5_4).' * F5;
    F3 = G3 * dV3 - ad(V3).' * G3 * V3 + computeAdjoint(T4_3).' * F4;
    F2 = G2 * dV2 - ad(V2).' * G2 * V2 + computeAdjoint(T3_2).' * F3;
    F1 = G1 * dV1 - ad(V1).' * G1 * V1 + computeAdjoint(T2_1).' * F2;

    % Torques
    torque1 = F1.' * A1;
    torque2 = F2.' * A2;
    torque3 = F3.' * A3;
    torque4 = F4.' * A4;
    
    %% Closed-form
    V_all = [V1;V2;V3;V4];

    F_all = [F1;F2;F3;F4];

    ad_V_all = blkdiag(ad(V1), ad(V2), ad(V3), ad(V4));
    
    ad_A_all_qd = blkdiag(ad(A1*q(1)), ad(A2*q(2)), ad(A3*q(3)), ad(A4*q(4)));
    
    % Build W(q)
    W = zeros(24, 24); 
    
    W(7:12, 1:6) = computeAdjoint(T2_1); 
    W(13:18, 7:12) = computeAdjoint(T3_2); 
    W(19:24, 13:18) = computeAdjoint(T4_3);
    
    % Build L(q)
    I = eye(24);
    L = inv(I - W);
    
    % Build Vbase,Vbase_dot
    Vbase = zeros(24,1);
    Vbase(1:6) = computeAdjoint(T1_0)*V0;

    Vbase_dot = zeros(24,1);
    Vbase_dot(1:6) = computeAdjoint(T1_0)*dV0;

    % Build Ftip
    F_tip = zeros(24,1);
    F_tip(end-5:end) = computeAdjoint(T5_4).' * F5;
    
    Vall = L * (A_all*qd + Vbase);
    Vall_dot = L * (A_all * qdd - ad_A_all_qd * W * Vall - ad_A_all_qd * Vbase + Vbase_dot);
    Fall = L' * (G_all * Vall_dot - ad_V_all' * G_all * Vall + F_tip);

    % Dynamic quantities
    torque = A_all' * Fall;
    M = A_all' * L' * G_all * L * A_all;
    C = -A_all' * L' * (G_all * L * ad_A_all_qd * W + ad_V_all' * G_all) * L * A_all * qd;
    g = A_all' * L' * G_all * L * Vbase_dot;
    
    % Forward Dynamics
    qdd = M \ (torque - C - g);

    % Euler Integration
    qd = qd + qdd * dt;
    q = q + qd * dt;

    % Visualization
    show(robot, q, 'PreservePlot', false);
    drawnow;
end

% Helper function to compute adjoint matrix
function Ad_T = computeAdjoint(T)
    R = T(1:3, 1:3); % Rotation
    p = T(1:3, 4);   % Translation
    Ad_T = [R, zeros(3); skew(p) * R, R];
end

% Helper function for skew-symmetric matrix
function S = skew(v)
    S = [  0,   -v(3),  v(2);
          v(3),   0,   -v(1);
         -v(2),  v(1),   0];
end

% Helper function to compute adjoint operator
function rad = ad(V)
    omega = V(1:3); % Angular velocity
    v = V(4:6);     % Linear velocity
    rad = [skew(omega), zeros(3); skew(v), skew(omega)];
end

% Convert a 6x1 screw axis vector into a 4x4 matrix representation
function S_hat = screwToMatrix(S)  
    % Extract rotational (omega) and translational (v) components
    omega = S(1:3); % First 3 components (rotation axis)
    v = S(4:6);     % Last 3 components (translation vector)
    
    % Compute the skew-symmetric matrix of omega
    omega_hat = [  0       -omega(3)  omega(2);
                   omega(3)  0       -omega(1);
                  -omega(2)  omega(1)  0      ];
    
    % Construct the 4x4 matrix
    S_hat = [omega_hat, v; 
             0, 0, 0, 0];
end