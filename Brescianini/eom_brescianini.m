function Xdot = eom_brescianini(t, X, k, param)

e3 = [0, 0, 1]';
m = param.m;
J = param.J;

[~, v, R, W] = split_to_states(X);

desired = command(t);
[f, M, ~, ~] = position_control(X, desired, k, param);

[f, M] = saturate_fM(f, M, param);

x_dot = v;
v_dot = param.g * e3 ...
    - f / m * R * e3;
if param.use_disturbances
     v_dot = v_dot + param.W_x * param.theta_x / m;
end

if param.use_disturbances
    W_dot = J \ (-hat(W) * J * W + M + param.W_R * param.theta_R);
else
    W_dot = J \ (-hat(W) * J * W + M);
end
R_dot = R * hat(W);

Xdot = [x_dot; v_dot; W_dot; reshape(R_dot,9,1)];
end