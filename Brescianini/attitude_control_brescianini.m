function [tau, error] = attitude_control_brescianini(...
    R, w, ...  % states
    Rd, wd, wd_dot, ...  % desired values
    param ...  % gains and parameters
)
%% Unpack parameters
J = param.J;
kp_xy = param.kp_xy;
kp_z = param.kp_z;
kd_xy = param.kd_xy;
kd_z = param.kd_z;
%%
q = quaternion(R, 'rotmat', 'frame');
qd = quaternion(Rd, 'rotmat', 'frame');

qe = qd*conj(q);

wd_bar = rotmat(qe, 'frame')'*wd;
we = wd_bar - w;

wd_bar_dot = hat(we)*wd_bar + rotmat(qe, 'frame')'*wd_dot;

qe = conj(qe);
[q0, q1, q2, q3] = parts(qe);

q0q3 = sqrt(q0^2 + q3^2);
B = [q0^2 + q3^2;
    q0*q1 - q2*q3;
    q0*q2 + q1*q3;
    0];
qe_red = B / q0q3;
qe_yaw = [q0; 0; 0; q3] / q0q3;

tilde_qe_red = qe_red(2:4);
tilde_qe_yaw = qe_yaw(2:4);

tau_ff = J*wd_bar_dot - hat(J*w)*w;

Kd = diag([kd_xy, kd_xy, kd_z]);
tau = kp_xy*tilde_qe_red ...
    + kp_z*sign(q0)*tilde_qe_yaw ...
    + Kd*we ...
    + tau_ff;
%%
error.qe = qe;
error.we = we;
end