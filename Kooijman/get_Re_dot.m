function Re_dot = get_Re_dot(v, v_dot, v_bar, v_bar_dot)

den = s(v_bar, v);
den_dot = s_dot(v, v_bar, v_dot, v_bar_dot);

num = hat(v_bar)*v;
num_dot = hat(v_bar_dot)*v + hat(v_bar)*v_dot;
Rrd1 = diff_num_den(num, num_dot, den, den_dot);
if norm(den) < 1e-3
    Rrd1 = Rrd1*0;
end

num = -(eye(3) - v_bar*v_bar')*v;
num_dot = (v_bar_dot*v_bar' + v_bar*v_bar_dot')*v ...
    - (eye(3) - v_bar*v_bar')*v_dot;
Rrd2 = diff_num_den(num, num_dot, den, den_dot);
if norm(den) < 1e-3
    Rrd2 = Rrd2*0;
end

Rrd3 = v_bar_dot;

Re_dot = [Rrd1, Rrd2, Rrd3];

end