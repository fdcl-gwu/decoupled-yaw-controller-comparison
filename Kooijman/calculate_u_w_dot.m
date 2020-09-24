function u_w_dot = calculate_u_w_dot(r1, r1_dot, ...
    r2, r2_dot, ...
    r3, r3_dot, ...
    r1_bar, r1_bar_dot, r1_bar_2dot, ...
    r3_bar, r3_bar_dot, ...
    param)

k2 = param.k2;

[R_e, R_r] = get_Re_Rr(r3, r3_bar);
R_r_dot = get_Rr_dot(r3, r3_dot, r3_bar, r3_bar_dot);
R_e_dot = get_Re_dot(r3, r3_dot, r3_bar, r3_bar_dot);

w1 = r1'*R_r*R_e'*r1_bar;
w1_dot = r1_dot'*R_r*R_e'*r1_bar ...
    + r1'*R_r_dot*R_e'*r1_bar ...
    + r1'*R_r*R_e_dot'*r1_bar ...
    + r1'*R_r*R_e'*r1_bar_dot;

w2 = r2'*R_r*R_e'*r1_bar;
w2_dot = r2_dot'*R_r*R_e'*r1_bar ...
    + r2'*R_r_dot*R_e'*r1_bar ...
    + r2'*R_r*R_e_dot'*r1_bar ...
    + r2'*R_r*R_e'*r1_bar_dot;

if abs(w1) > abs(w2)
    theta2 = r2'*R_r*R_e'*r1_bar_dot;
    theta2_dot = r2_dot'*R_r*R_e'*r1_bar_dot ...
        + r2'*R_r_dot*R_e'*r1_bar_dot ...
        + r2'*R_r*R_e_dot'*r1_bar_dot ...
        + r2'*R_r*R_e'*r1_bar_2dot;
    
    num = theta2;
    num_dot = theta2_dot;

    den = w1;
    den_dot = w1_dot;

    u_w_FF_dot = diff_num_den(num, num_dot, den, den_dot);
else
    theta1 = -r1'*R_r*R_e'*r1_bar_dot;
    theta1_dot = -r1_dot'*R_r*R_e'*r1_bar_dot ...
        - r1'*R_r_dot*R_e'*r1_bar_dot ...
        - r1'*R_r*R_e_dot'*r1_bar_dot ...
        - r1'*R_r*R_e'*r1_bar_2dot;

    num = theta1;
    num_dot = theta1_dot;

    den = w2;
    den_dot = w2_dot;

    u_w_FF_dot = diff_num_den(num, num_dot, den, den_dot);
end

if w1 >= 0
    u_w_FB_dot = k2*w2_dot;
elseif w1 < 0 && w2 < 0
    u_w_FB_dot = -k2;
else
    u_w_FB_dot = k2;
end

u_w_dot = u_w_FB_dot + u_w_FF_dot;

end