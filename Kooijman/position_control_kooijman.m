function [T, tau, error, calculated] = position_control_kooijman( ...
    X, desired, param)

[x, v, R, W] = split_to_states(X);
R_dot = R*hat(W);

m = param.m;
J = param.J;
g = param.g;

kp = param.kp;
kv = param.kv;
kW = param.kW;
kwy = param.kwy;

e3 = [0, 0, 1]';

%%
error.x = desired.x - x;
error.v = desired.v - v;

%%
r1 = R(:,1);
r2 = R(:,2);
r3 = R(:,3);

r1_dot = R_dot(:,1);
r2_dot = R_dot(:,2);
r3_dot = R_dot(:,3);

b3 = R*e3;
b3_dot = R_dot*e3;
%%
b = -desired.x_2dot + g*e3;

T_bar = m*norm(b);
T_msqrt3 = T_bar / (sqrt(3)*m);
L_lower = [-T_msqrt3, -T_msqrt3, g - T_msqrt3]';
L_upper = [T_msqrt3, T_msqrt3, g + T_msqrt3]';

u_bar = desired.x_2dot;
% u = saturate_u(u_bar + kp*error.x + kv*error.v, L_lower, L_upper);
u = u_bar + kp*error.x + kv*error.v;

a_bar_ref = g*e3 - u;
n_a_bar_ref = norm(a_bar_ref);

T = m*n_a_bar_ref;
u_bar_dot = desired.x_3dot;
v_dot = g * e3 ...
    - T / m * b3;
error.a = desired.x_2dot - v_dot;
%u_dot = saturate_u(u_bar_dot + kp*error.v + kv*error.a, L_lower, L_upper);
u_dot = u_bar_dot + kp*error.v + kv*error.a;
a_ref_dot = -u_dot;

n_a_ref_dot = a_bar_ref'*a_ref_dot / n_a_bar_ref;
T_dot = m*n_a_ref_dot;
v_2dot = - T_dot / m * b3 ...
    - T / m * b3_dot;
error.a_dot = desired.x_3dot - v_2dot;

u_bar_2dot = desired.x_4dot;
% u_2dot = saturate_u(u_bar_2dot + kp*error.a + kv*error.a_dot, ...
%     L_lower, L_upper);
u_2dot = u_bar_2dot + kp*error.a + kv*error.a_dot;
a_ref_2dot = -u_2dot;

[r3_bar, r3_bar_dot, r3_bar_2dot] = deriv_unit_vector(a_bar_ref, ...
    a_ref_dot, a_ref_2dot);

% phi_bar = atan2(desired.b1(2), desired.b1(1));
phi_bar = desired.yaw;
phi_bar_dot = desired.w;
phi_bar_2dot = desired.w_dot;

r_yaw = [-sin(phi_bar), cos(phi_bar), 0]';
r_yaw_dot = [-cos(phi_bar)*phi_bar_dot;
    -sin(phi_bar)*phi_bar_dot;
    0];
r_yaw_2dot = [sin(phi_bar)*phi_bar_dot^2 + -cos(phi_bar)*phi_bar_2dot;
    -cos(phi_bar)*phi_bar_dot^2 - sin(phi_bar)*phi_bar_2dot;
    0];

num = hat(r_yaw)*r3_bar;
num_dot = hat(r_yaw_dot)*r3_bar + hat(r_yaw)*r3_bar_dot;
num_2dot = hat(r_yaw_2dot)*r3_bar ...
    + hat(r_yaw_dot)*r3_bar_dot ...
    + hat(r_yaw_dot)*r3_bar_dot ...
    + hat(r_yaw)*r3_bar_2dot;

den = s(r_yaw, r3_bar);
den_dot = s_dot(r_yaw, r3_bar, r_yaw_dot, r3_bar_dot);
den_2dot = s_2dot(r_yaw, r3_bar, ...
    r_yaw_dot, r3_bar_dot, ...
    r_yaw_2dot, r3_bar_2dot);

r1_bar = num/den;
r1_bar_dot = diff_num_den(num, num_dot, den, den_dot);
r1_bar_2dot = diff2_num_den(num, num_dot, num_2dot, ...
    den, den_dot, den_2dot);

r2_bar = hat(r3_bar)*r1_bar;

u_v = calculate_u_v(r3, r3_bar, r3_bar_dot, r1, param);
u_w = calculate_u_w(r1, r2, r3, r1_bar, r1_bar_dot, r3_bar, param);

[R_e, R_r] = get_Re_Rr(r3, r3_bar);

% r3_dot = (eye(3) - r3*r3')*u_v;
R_r_dot = get_Rr_dot(r3, r3_dot, r3_bar, r3_bar_dot);
w_r = vee(R_r'*R_r_dot);

R_e_dot = get_Re_dot(r3, r3_dot, r3_bar, r3_bar_dot);
w_e = vee(R_e'*R_e_dot);

W_bar1 = -r2'*u_v;
W_bar2 = r1'*u_v;

if abs(r3'*r3_bar) > 1e-3
    w1 = r1'*R_r*R_e'*r1_bar;
    w2 = r2'*R_r*R_e'*r1_bar;
else
    w1 = r1'*r1_bar;
    w2 = r2'*r2_bar;
end

beta1 = w2*r3'*R_r*R_e'*r1_bar - r1'*R_r*hat(w_r - w_e)*R_e'*r1_bar;
beta2 = w1*r3'*R_r*R_e'*r1_bar + r2'*R_r*hat(w_r - w_e)*R_e'*r1_bar;

if abs(w1) > abs(w2)
    w_r = beta2/w1;
else
    w_r = beta1/w2;
end

W_bar = [W_bar1, W_bar2, u_w + w_r]';
%%
r3_dot = R_dot(:,3);
u_v_dot = calculate_u_v_dot(r3, r3_dot, ...
    r3_bar, r3_bar_dot, r3_bar_2dot, ...
    r1_dot, param);

u_w_dot = calculate_u_w_dot(r1, r1_dot, ...
    r2, r2_dot, ...
    r3, r3_dot, ...
    r1_bar, r1_bar_dot, r1_bar_2dot, ...
    r3_bar, r3_bar_dot, ...
    param);
%%
r3_2dot = (- r3_dot*r3' - r3*r3_dot')*u_v ...
    + (eye(3) - r3*r3')*u_v_dot;
%%
w1_dot = r1_dot'*R_r*R_e'*r1_bar ...
    + r1'*R_r_dot*R_e'*r1_bar ...
    + r1'*R_r*R_e_dot'*r1_bar ...
    + r1'*R_r*R_e'*r1_bar_dot;

w2_dot = r2_dot'*R_r*R_e'*r1_bar ...
    + r2'*R_r_dot*R_e'*r1_bar ...
    + r2'*R_r*R_e_dot'*r1_bar ...
    + r2'*R_r*R_e'*r1_bar_dot;

R_r_2dot = get_Rr_2dot(r3, r3_dot, r3_2dot, ...
    r3_bar, r3_bar_dot, r3_bar_2dot);
R_e_2dot = get_Re_2dot(r3, r3_dot, r3_2dot, ...
    r3_bar, r3_bar_dot, r3_bar_2dot);

w_r_dot = vee(R_r_dot'*R_r_dot) + vee(R_r'*R_r_2dot);
w_e_dot = vee(R_e_dot'*R_e_dot) + vee(R_e'*R_e_2dot);

beta1_dot = w2_dot*r3'*R_r*R_e'*r1_bar ...
    + w2*r3_dot'*R_r*R_e'*r1_bar ...
    + w2*r3'*R_r_dot*R_e'*r1_bar ...
    + w2*r3'*R_r*R_e_dot'*r1_bar ...
    + w2*r3'*R_r*R_e'*r1_bar_dot ...
    - r1_dot'*R_r*hat(w_r - w_e)*R_e'*r1_bar ...
    - r1'*R_r_dot*hat(w_r - w_e)*R_e'*r1_bar ...
    - r1'*R_r*hat(w_r_dot - w_e_dot)*R_e'*r1_bar ...
    - r1'*R_r*hat(w_r - w_e)*R_e_dot'*r1_bar ...
    - r1'*R_r*hat(w_r - w_e)*R_e'*r1_bar_dot;

beta2_dot = w1_dot*r3'*R_r*R_e'*r1_bar ...
    + w1*r3_dot'*R_r*R_e'*r1_bar ...
    + w1*r3'*R_r_dot*R_e'*r1_bar ...
    + w1*r3'*R_r*R_e_dot'*r1_bar ...
    + w1*r3'*R_r*R_e'*r1_bar_dot ...
    + r2_dot'*R_r*hat(w_r - w_e)*R_e'*r1_bar ...
    + r2'*R_r_dot*hat(w_r - w_e)*R_e'*r1_bar ...
    + r2'*R_r*hat(w_r_dot - w_e_dot)*R_e'*r1_bar ...
    + r2'*R_r*hat(w_r - w_e)*R_e_dot'*r1_bar ...
    + r2'*R_r*hat(w_r - w_e)*R_e'*r1_bar_dot;

if abs(w1) > abs(w2)
    w_r_dot = diff_num_den(beta2, beta2_dot, w1, w1_dot);
else
    w_r_dot = diff_num_den(beta1, beta1_dot, w2, w2_dot);
end

W1_dot = - r2_dot'*u_v - r2'*u_v_dot;
W2_dot = r1_dot'*u_v + r1'*u_v_dot;
W3_dot = u_w_dot + w_r_dot;
W_bar_dot = [W1_dot, W2_dot, W3_dot]';
%%
Rd = [r1_bar, r2_bar, r3_bar];

Wd = W_bar;
Wd_dot = W_bar_dot;

kW = diag([kW, kW, kwy]);

eW = W - Wd;
tau = -kW*eW + hat(W)*J*W + J*Wd_dot;

%% Saving data
calculated.b3 = r3_bar;
calculated.b3_dot = r3_bar_dot;
calculated.R = Rd;
error.R = 0.5*vee(Rd'*R - R'*Rd);
end