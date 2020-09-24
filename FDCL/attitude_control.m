function [M, dot_theta_R, eb, ew, ey, ewy] = attitude_control(...
    R, W, bar_theta_R, ...  % states
    b3d, b3d_dot, b3d_ddot, b1c, wc3, wc3_dot, ...  % desired values
    k, param ...  % gains and parameters
)
%% Unpack other parameters
J = param.J;
c2 = param.c2;
c3 = param.c3;

gamma_R = param.gamma_R;
W_R = param.W_R;

W_R_1 = W_R(1,:);
W_R_2 = W_R(2,:);
W_R_3 = W_R(3,:);
%% Body axes
e1 = [1, 0, 0]';
e2 = [0, 1 ,0]';
e3 = [0, 0, 1]';

b1 = R * e1;
b2 = R * e2;
b3 = R * e3;
%% Roll/pitch dynamics
kb = k.R;
kw = k.W;

w = W(1) * b1 + W(2) * b2;
b3_dot = hat(w) * b3;

wd = hat(b3d) * b3d_dot;
wd_dot = hat(b3d) * b3d_ddot;

eb = hat(b3d) * b3;
ew = w + hat(b3)^2 * wd;
tau = - kb * eb ...
    - kw * ew ...
    - J(1,1) * dot(b3, wd) * b3_dot ...
    - J(1,1) * hat(b3)^2 * wd_dot ...
    - W_R_1 * bar_theta_R * b1 - W_R_2 * bar_theta_R * b2;

tau1 = dot(b1, tau);
tau2 = dot(b2, tau);

M1 = tau1 + J(3,3) * W(3) * W(2);
M2 = tau2 - J(3,3) * W(3) * W(1);
%% Yaw dynamics
ey = -dot(b2, b1c);
ewy = W(3) - wc3;

M3 = - k.y * ey ...
    - k.wy * ewy ...
    - W_R_3 * bar_theta_R ...
    + J(3,3) * wc3_dot;
%% Attitude adaptive term
ew_c2eb = ew + c2 * eb;
dot_theta_R = gamma_R * W_R_1' * ew_c2eb' * b1 ...
    + gamma_R * W_R_2' * ew_c2eb' * b2 ...
    + gamma_R * W_R_3' * (ewy + c3 * ey);

M = [M1, M2, M3]';
end