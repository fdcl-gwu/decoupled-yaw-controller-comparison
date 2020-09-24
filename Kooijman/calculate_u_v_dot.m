function u_v_dot = calculate_u_v_dot(v, v_dot, ...
    v_bar, v_bar_dot, v_bar_2dot, ...
    r1_dot, param)

k1 = param.k1;

if (v'*v_bar >= 0)
    u_v_FB_dot = k1*v_bar_dot;
elseif (v == -v_bar)
    u_v_FB_dot = k1*r1_dot;
else
    num = k1*v_bar;
    num_dot = k1*v_bar_dot;
    den = s(v, v_bar);
    den_dot = s_dot(v, v_bar, v_dot, v_bar_dot);
    u_v_FB_dot = diff_num_den(num, num_dot, den, den_dot);
end

if (v == v_bar)
    u_v_FF_dot = v_bar_2dot;
elseif (v == -v_bar)
    u_v_FF_dot = -v_bar_2dot;
else
    vxvbar = cross(v, v_bar);
    vxvbar_dot = cross(v_dot, v_bar) + cross(v, v_bar_dot);

    num = vxvbar*vxvbar' - (eye(3) - v*v')*v_bar*v';
    num_dot = vxvbar_dot*vxvbar' ...
        + vxvbar*vxvbar_dot' ...
        - (- v_dot*v' - v*v_dot')*v_bar*v' ...
        - (eye(3) - v*v')*v_bar_dot*v' ...
        - (eye(3) - v*v')*v_bar*v_dot';

    den = s(v, v_bar)^2;
    den_dot = 2*s(v, v_bar)*s_dot(v, v_bar, v_dot, v_bar_dot);

    theta = num / den;
    theta_dot = diff_num_den(num, num_dot, den, den_dot);

    u_v_FF_dot = theta_dot*v_bar_dot + theta*v_bar_2dot;
end

u_v_dot = u_v_FB_dot + u_v_FF_dot;

end