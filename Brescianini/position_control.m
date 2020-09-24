function [f, M, error, calculated] ...
    = position_control(X, desired, k, param)
% Unpack states
[x, v, R, W] = split_to_states(X);

m = param.m;
g = param.g;
e3 = [0, 0, 1]';

error.x = x - desired.x;
error.v = v - desired.v;
A = - k.x * error.x ...
    - k.v * error.v ...
    - m * g * e3 ...
    + m * desired.x_2dot;

b3 = R * e3;
f = -dot(A, b3);
ea = g * e3 ...
    - f / m * b3 ...
    - desired.x_2dot;
A_dot = - k.x * error.v ...
    - k.v * ea ...
    + m * desired.x_3dot;

b3_dot = R * hat(W) * e3;
f_dot = -dot(A_dot, b3) - dot(A, b3_dot);
eb = - f_dot / m * b3 - f / m * b3_dot - desired.x_3dot;
A_ddot = - k.x * ea ...
    - k.v * eb ...
    + m * desired.x_4dot;

[b3c, b3c_dot, b3c_ddot] = deriv_unit_vector(-A, -A_dot, -A_ddot);

A2 = -hat(desired.b1) * b3c;
A2_dot = -hat(desired.b1_dot) * b3c - hat(desired.b1) * b3c_dot;
A2_ddot = - hat(desired.b1_2dot) * b3c ...
    - 2 * hat(desired.b1_dot) * b3c_dot ...
    - hat(desired.b1) * b3c_ddot;

[b2c, b2c_dot, b2c_ddot] = deriv_unit_vector(A2, A2_dot, A2_ddot);

b1c = hat(b2c) * b3c;
b1c_dot = hat(b2c_dot) * b3c+hat(b2c)*b3c_dot;
b1c_ddot = hat(b2c_ddot) * b3c ...
    + 2 * hat(b2c_dot) * b3c_dot ...
    + hat(b2c) * b2c_ddot;

Rc = [b1c, b2c, b3c];
Rc_dot = [b1c_dot, b2c_dot, b3c_dot];
Rc_ddot = [b1c_ddot, b2c_ddot, b3c_ddot];

Wc = vee(Rc' * Rc_dot);
Wc_dot = vee(Rc' * Rc_ddot - hat(Wc)^2);

%% Run attitude controller
[M, ~] = attitude_control_brescianini(...
    R, W, ...  % states
    Rc, Wc, Wc_dot, ...  % desired values
    param ...  % gains and parameters
);
error.R = 0.5*vee(Rc'*R - R'*Rc);
%% Saving data
calculated.b3 = b3c;
calculated.b3_dot = b3c_dot;
calculated.b3_ddot = b3c_ddot;
calculated.b1 = b1c;
calculated.R = Rc;
calculated.W = Wc;
calculated.W_dot = Wc_dot;
end