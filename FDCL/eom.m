function Xdot = eom(t, X, k, param)

e3 = [0, 0, 1]';
m = param.m;
J = param.J;

W_x = param.W_x;
W_R = param.W_R;

[~, v, R, W, ~, ~] = split_to_states(X);

desired = command(t);
[f, M, bar_theta_x_dot, bar_theta_R_dot, ~, ~] = position_control(X, ...
    desired, k, param);

[f, M] = saturate_fM(f, M, param);

xdot = v;
vdot = param.g * e3 ...
    - f / m * R * e3 + W_x * param.theta_x / m;
Wdot = J \ (-hat(W) * J * W + M + W_R * param.theta_R);
Rdot = R * hat(W);

if ~param.use_disturbances
    bar_theta_x_dot = 0*bar_theta_x_dot;
    bar_theta_R_dot = 0*bar_theta_R_dot;
end

Xdot=[xdot; vdot; Wdot; reshape(Rdot,9,1); ...
    bar_theta_x_dot; bar_theta_R_dot];
end