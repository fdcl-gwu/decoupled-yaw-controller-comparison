function [f, M, bar_theta_x_dot, theta_R_dot, error, calculated] = ...
    position_control(X, desired, k, param)

[x, v, R, W, bar_theta_x, bar_theta_R] = split_to_states(X);

c1 = param.c1;
m = param.m;
g = param.g;

W_x = param.W_x;
W_x_dot = zeros(3);
W_x_2dot = zeros(3);

e3 = [0, 0, 1]';

%%
error.x = x - desired.x;
error.v = v - desired.v;
A = - k.x * error.x ...
    - k.v * error.v ...
    - m * g * e3 ...
    + m * desired.x_2dot ...
    - W_x * bar_theta_x;
%%
gamma_x = param.gamma_x;
c1 = param.c1;
ev_c1ex = error.v + c1 * error.x;

norm_theta_x = norm(bar_theta_x);
if norm_theta_x < param.B_theta_x || ...
    (norm_theta_x == param.B_theta_x && bar_theta_x' * W_x' * ev_c1ex <= 0)
    
    bar_theta_x_dot = gamma_x * W_x' * ev_c1ex;
else
    I_theta = eye(3) ...
        - bar_theta_x * bar_theta_x' / (bar_theta_x' * bar_theta_x);
    bar_theta_x_dot = gamma_x * I_theta * W_x' * ev_c1ex;
end
%%
b3 = R * e3;
f = -dot(A, b3);
ev_dot = g * e3 ...
    - f / m * b3 ...
    - desired.x_2dot ...
    + W_x * bar_theta_x / m;
A_dot = - k.x * error.v ...
    - k.v * ev_dot ...
    + m * desired.x_3dot ...
    - W_x_dot * bar_theta_x ...
    - W_x * bar_theta_x_dot;
%%
norm_theta_x = norm(bar_theta_x);
if norm_theta_x < param.B_theta_x || ...
    (norm_theta_x == param.B_theta_x && bar_theta_x' * W_x' * ev_c1ex <= 0)
    
    bar_theta_x_2dot = gamma_x * W_x_dot' * ev_c1ex ...
        + gamma_x * W_x' * (ev_dot + c1 * error.v);
else
    I_theta = eye(3) ...
        - bar_theta_x * bar_theta_x' / (bar_theta_x' * bar_theta_x);
    
    num = norm_theta_x * (bar_theta_x_dot * bar_theta_x' ...
        + bar_theta_x * bar_theta_x_dot') ...
        - 2 * (bar_theta_x * bar_theta_x') * bar_theta_x_dot;
    I_theta_dot = - num / norm_theta_x^3;
    bar_theta_x_2dot = gamma_x * I_theta_dot * W_x' * ev_c1ex ...
        + gamma_x * I_theta * W_x_dot' * ev_c1ex ...
        + gamma_x * I_theta * W_x' * (ev_dot + c1 * error.v);
end
%%
b3_dot = R * hat(W) * e3;
f_dot = -dot(A_dot, b3) - dot(A, b3_dot);
ev_2dot = - f_dot / m * b3 - f / m * b3_dot - desired.x_3dot ...
    + W_x_dot * bar_theta_x / m + W_x * bar_theta_x_dot / m;
A_ddot = - k.x * ev_dot ...
    - k.v * ev_2dot ...
    + m * desired.x_4dot ...
    - W_x_2dot * bar_theta_x ...
    - 2 * W_x_dot * bar_theta_x_dot ...
    - W_x * bar_theta_x_2dot;
%%
[b3c, b3c_dot, b3c_ddot] = deriv_unit_vector(-A, -A_dot, -A_ddot);

A2 = -hat(desired.b1) * b3c;
A2_dot = -hat(desired.b1_dot) * b3c - hat(desired.b1) * b3c_dot;
A2_ddot = - hat(desired.b1_2dot) * b3c ...
    - 2 * hat(desired.b1_dot) * b3c_dot ...
    - hat(desired.b1) * b3c_ddot;

[b2c, b2c_dot, b2c_ddot] = deriv_unit_vector(A2, A2_dot, A2_ddot);

b1c = hat(b2c) * b3c;
b1c_dot = hat(b2c_dot) * b3c + hat(b2c)*b3c_dot;
b1c_ddot = hat(b2c_ddot) * b3c ...
    + 2 * hat(b2c_dot) * b3c_dot ...
    + hat(b2c) * b2c_ddot;
%%
Rc = [b1c, b2c, b3c];
Rc_dot = [b1c_dot, b2c_dot, b3c_dot];
Rc_ddot = [b1c_ddot, b2c_ddot, b3c_ddot];
%%
Wc = vee(Rc' * Rc_dot);
Wc_dot = vee(Rc' * Rc_ddot - hat(Wc)^2);
%%
W3 = dot(R * e3, Rc * Wc);
W3_dot = dot(R * e3, Rc * Wc_dot) ...
    + dot(R * hat(W) * e3, Rc * Wc);
%% Run attitude controller
if param.use_decoupled_controller
    [M, theta_R_dot, error.R, error.W, error.y, error.Wy] ...
        = attitude_control( ...
        R, W, bar_theta_R, ...
        b3c, b3c_dot, b3c_ddot, b1c, W3, W3_dot, ...
        k, param);
    
    % For comparison with non-decoupled controller
    error.R = 1 / 2 * vee(Rc' * R - R' * Rc);
else
    [M, theta_R_dot, error.R, error.W] = attitude_control_coupled(...
        R, W, bar_theta_R, Rc, Wc, Wc_dot, k, param);
end
%% Saving data
calculated.b3 = b3c;
calculated.b3_dot = b3c_dot;
calculated.b3_ddot = b3c_ddot;
calculated.b1 = b1c;
calculated.R = Rc;
calculated.W = Wc;
calculated.W_dot = Wc_dot;
calculated.W3 = dot(R * e3, Rc * Wc);
calculated.W3_dot = dot(R * e3, Rc * Wc_dot) ...
    + dot(R * hat(W) * e3, Rc * Wc);
end