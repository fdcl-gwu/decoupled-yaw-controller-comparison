function Xdot = eom_kooijman(t, X, param)

e3 = [0, 0, 1]';
m = param.m;
J = param.J;

[~, v, R, W] = split_to_states(X);

desired = command(t);
[T, tau, ~, ~] = position_control_kooijman(X, desired, param);

thr = fM_to_thr(T, tau, param)';
[T, tau] = saturate_fM(T, tau, param);

x_dot = v;
v_dot = param.g*e3 - T*R*e3 / m;
if param.use_disturbances
    v_dot = v_dot + param.W_x*param.theta_x/m;
end

R_dot = R * hat(W);

if param.use_disturbances
    W_dot = J \ (-hat(J*W)*W + tau + param.W_R*param.theta_R);
else
    W_dot = J \ (-hat(W)*J*W + tau);
end

Xdot = [x_dot; v_dot; reshape(R_dot,9,1); W_dot];
end