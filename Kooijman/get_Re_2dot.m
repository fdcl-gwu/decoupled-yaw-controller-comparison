function Re_2dot = get_Re_2dot(v, v_dot, v_2dot, ...
    v_bar, v_bar_dot, v_bar_2dot)

den = s(v_bar, v);
den_dot = s_dot(v, v_bar, v_dot, v_bar_dot);
den_2dot = s_2dot(v, v_bar, v_dot, v_bar_dot, v_2dot, v_bar_2dot);
%%
num = hat(v_bar)*v;
num_dot = hat(v_bar_dot)*v + hat(v_bar)*v_dot;
num_2dot = hat(v_bar_2dot)*v + hat(v_bar_dot)*v_dot ...
    + hat(v_bar_dot)*v_dot + hat(v_bar)*v_2dot;
Rrd1 = diff2_num_den(num, num_dot, num_2dot, ...
    den, den_dot, den_2dot);
%%
num1 = (eye(3) - v_bar*v_bar');
num1_dot = -(v_bar_dot*v_bar' + v_bar*v_bar_dot');
num1_2dot = -(v_bar_2dot*v_bar' + v_bar_dot*v_bar_dot' ...
    + v_bar_dot*v_bar_dot' + v_bar*v_bar_2dot');

num = -num1*v;
num_dot = -num1_dot*v - num1*v_dot;
num_2dot = - num1_2dot*v - 2*(num1_dot*v_dot) - num1*v_2dot;

Rrd2 = diff2_num_den(num, num_dot, num_2dot, ...
    den, den_dot, den_2dot);
%%
Rrd3 = v_bar_2dot;
%%
Re_2dot = [Rrd1, Rrd2, Rrd3];

end