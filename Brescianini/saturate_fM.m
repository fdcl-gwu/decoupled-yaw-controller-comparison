function [f_hat, tau_hat] = saturate_fM(f, tau, param)

l = param.d;
c_tf = param.c_tf;

u_max = 8;
u_min = 0.1;

tau_hat = [0, 0, 0]';
tau_max_xy = (u_max - u_min)*l;
for i = 1:2
    tau_hat(i) = saturate(tau(i), -tau_max_xy, tau_max_xy);
end

tau_hat_x = tau_hat(1);
tau_hat_y = tau_hat(2);
f_min = 4*u_min + abs(tau_hat_x)/l + abs(tau_hat_y)/l;
f_max = 4*u_max - abs(tau_hat_x)/l - abs(tau_hat_y)/l;
f_hat = saturate(f, f_min, f_max);

tau_min_z_list = c_tf*[4*u_min - f_hat + 2*abs(tau_hat_x)/l;
    -4*u_max + f_hat + 2*abs(tau_hat_y)/l];
tau_min_z = max(tau_min_z_list);

tau_max_z_list = c_tf*[4*u_max - f_hat - 2*abs(tau_hat_x)/l;
    -4*u_min + f_hat - 2*abs(tau_hat_y)/l];
tau_max_z = min(tau_max_z_list);

tau_hat(3) = saturate(tau(3), tau_min_z, tau_max_z);
end