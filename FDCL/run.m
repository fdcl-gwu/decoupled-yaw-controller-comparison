close all;
clear all;
addpath('../Common');
%% Simulation mode
% Uncomment to use geometric adaptive decoupled-yaw controller.
param.use_decoupled_controller = true;

% Uncomment to use geometric adaptive coupled-yaw controller in
% Geometric adaptive tracking control of a quadrotor unmanned aerial 
% vehicle on {SE(3)}, F. Goodarzi, D. Lee, and T. Lee.
% param.use_decoupled_controller = false;
%% Disturbances
% Uncomment to use disturbances.
% param.use_disturbances = true;

% Uncomment to remove disturbances.
param.use_disturbances = false;
%% Simulation parameters
t = 0:0.01:10;
N = length(t);

% Quadrotor
J1 = 0.02;
J2 = 0.02;
J3 = 0.04;
param.J = diag([J1, J2, J3]);
param.m = 2;
param.g = 9.81;

param.d = 0.169;
param.c_tf = 0.0135;

% Fixed disturbance
if param.use_disturbances
    param.W_x = eye(3);
    param.theta_x = [1, 0.8, -1]';

    param.W_R = eye(3);
    param.theta_R = [0.1, 0.1, -0.1]';
else
    param.W_x = zeros(3);
    param.theta_x = [0, 0, 0]';

    param.W_R = eye(3);
    param.theta_R = [0, 0, 0]';
end
%% Controller gains
k.x = 12;
k.v = 8;
k.R = 6;
k.W = 2;
k.y = 2;
k.wy = 0.8;
param.gamma_x = 2;
param.gamma_R = 10;
param.B_theta_x = 10;
%%
param.c1 = min(sqrt(k.x / param.m), ...
    4 * k.x * k.v / (k.v ^2 + 4 * param.m * k.x));

B2 = 1;
J1 = param.J(1,1);
param.c2 = min(sqrt(k.R / J1), ...
    4 * k.R * k.W / ((k.W + J1 * B2)^2 + 4 * J1 * k.R));

J3 = param.J(3,3);
param.c3 = min(sqrt(k.y / J3), ...
    4 * k.y * k.wy / (k.wy^2 + 4 * J1 * k.R));
%% Initial conditions
x0 = [1, -1, 0]'; % for circle trajectory
% x0 = [0, 0, 0]'; % for line trajectory

v0 = [0, 0, 0]';
e3 = [0, 0, 1]';
R0 = expm((pi - 0.01) * hat(e3));
W0 = [0, 0, 0]';

X0 = [x0; v0; W0; reshape(R0,9,1); zeros(6,1)];
%% Numerical integration
[t, X] = ode45(@(t, XR) eom(t, XR, k, param), t, X0, ...
    odeset('RelTol', 1e-6, 'AbsTol', 1e-6));
%% Output arrays
% Create empty arrays to save data
[e, d, R, f, M] = generate_output_arrays(N);
%% Post processing
x = X(:, 1:3)';
v = X(:, 4:6)';
W = X(:, 7:9)';
theta_x = X(:, 19:21)';
theta_R = X(:, 22:24)';

b1 = zeros(3, N);
b1c = zeros(3, N);

thr = zeros(4, N);

avg_ex = 0;
avg_eR = 0;
avg_f = 0;

converge_t = 0;
is_converged = false;
converge_ex = 0.02;

for i = 1:N
    R(:,:,i) = reshape(X(i,10:18), 3, 3);
    
    des = command(t(i));
    [f(i), M(:,i), ~, ~, err, calc] = position_control(X(i,:)', des, ...
        k, param);
    
    % Unpack errors
    e.x(:,i) = err.x;
    e.v(:,i) = err.v;
    e.R(:,i) = err.R;
    e.W(:,i) = err.W;
    
    if param.use_decoupled_controller
        e.y(i) = err.y;
        e.Wy(i) = err.Wy;
    end
    
    [f(i), M(:,i)] = saturate_fM(f(i), M(:,i), param);
    thr(:,i) = fM_to_thr(f(i), M(:,i), param);
    
    % Unpack desired values
    d.x(:,i) = des.x;
    d.v(:,i) = des.v;
    d.b1(:,i) = des.b1;
    d.R(:,:,i) = calc.R;
    b1(:,i) = R(:,:,i) * [1, 0, 0]';
    b1c(:,i) = calc.b1;
    
    norm_ex = norm(err.x);
    norm_eR = norm(err.R);
    
    % Find normalized errors
    avg_ex = avg_ex + norm_ex;
    avg_eR = avg_eR + norm_eR;
    
    norm_f = norm(thr(:,i));
    avg_f = avg_f + norm_f;
    
    if norm_ex < converge_ex
        if ~is_converged
            converge_t = t(i);
            is_converged = true;
        end
    end
end

avg_ex = avg_ex / N
avg_eR = avg_eR / N
avg_f = avg_f / N
converge_t

%% Plots
linetype = 'k';
linewidth = 1;
xlabel_ = 'time (s)';

figure(1);
plot_3x1(t, e.R, '', xlabel_, 'e_R', linetype, linewidth)
set(gca, 'FontName', 'Times New Roman');

figure(2);
plot_3x1(t, e.x, '', xlabel_, 'e_x', linetype, linewidth)
set(gca, 'FontName', 'Times New Roman');

figure(3);
plot_3x1(t, x, '', xlabel_, 'x', linetype, linewidth)
plot_3x1(t, d.x, '', xlabel_, 'x', 'r:', linewidth)
set(gca, 'FontName', 'Times New Roman');

figure(4);
plot_3x1(t, theta_x - param.theta_x, '', xlabel_, ...
    '\tilde\theta_x', linetype, linewidth)
set(gca, 'FontName', 'Times New Roman');

figure(5);
plot_3x1(t, theta_R - param.theta_R, '', xlabel_, ...
    '\tilde\theta_R', linetype, linewidth)
set(gca, 'FontName', 'Times New Roman');

figure(6);
plot_4x1(t, thr, '', xlabel_, 'f', linetype, linewidth)
set(gca, 'FontName', 'Times New Roman');

%% Save data
if param.use_decoupled_controller
    save('decoupled.mat');
else
    save('coupled.mat');
end