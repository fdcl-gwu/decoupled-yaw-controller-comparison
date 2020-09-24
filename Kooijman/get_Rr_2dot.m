function Rr_2dot = get_Rr_2dot(v, v_dot, v_2dot, ...
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
num = (eye(3) - v*v')*v_bar;

num1 = v_dot*v' + v*v_dot';
num1_dot = v_2dot*v' + 2*(v_dot*v_dot') + v*v_2dot';

num_dot = -num1*v_bar + (eye(3) - v*v')*v_bar_dot;
num_2dot = -num1_dot*v_bar ...
    - num1*v_bar_dot ...
    + (- v_dot*v' - v*v_dot')*v_bar_dot ...
    + (eye(3) - v*v')*v_bar_2dot;
Rrd2 = diff2_num_den(num, num_dot, num_2dot, ...
    den, den_dot, den_2dot);
%%
Rrd3 = v_2dot;
%%
Rr_2dot = [Rrd1, Rrd2, Rrd3];

end