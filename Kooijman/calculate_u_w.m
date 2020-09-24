function u_w = calculate_u_w(r1, r2, r3, r1_bar, r1_bar_dot, ...
    r3_bar, param)

k2 = param.k2;

[R_e, R_r] = get_Re_Rr(r3, r3_bar);

w1 = r1'*R_r*R_e'*r1_bar;
w2 = r2'*R_r*R_e'*r1_bar;

if abs(w1) > abs(w2)
    theta2 = r2'*R_r*R_e'*r1_bar_dot;
    u_w_FF = theta2 / w1;
else
    theta1 = -r1'*R_r*R_e'*r1_bar_dot;
    u_w_FF = theta1 / w2;
end

if w1 >= 0
    u_w_FB = k2*w2;
elseif w1 < 0 && w2 < 0
    u_w_FB = -k2;
else
    u_w_FB = k2;
end

u_w = u_w_FB + u_w_FF;
end