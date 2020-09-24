function u_v = calculate_u_v(v, v_bar, v_bar_dot, r1, param)

k1 = param.k1;

if (v'*v_bar >= 0)
    u_v_FB = k1*v_bar;
elseif (v == -v_bar)
    u_v_FB = k1*r1;
else
    u_v_FB = k1*v_bar / s(v, v_bar);
end

if (v == v_bar)
    u_v_FF = v_bar_dot;
elseif (v == -v_bar)
    u_v_FF = -v_bar_dot;
else
    vxvbar = cross(v, v_bar);
    theta = 1 / s(v, v_bar)^2 ...
        * (vxvbar*vxvbar' - (eye(3) - v*v')*v_bar*v');
    u_v_FF = theta*v_bar_dot;
end

u_v = u_v_FB + u_v_FF;

end