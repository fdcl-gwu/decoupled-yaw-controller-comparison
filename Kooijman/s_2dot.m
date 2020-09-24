function x = s_2dot(a, b, a_dot, b_dot, a_2dot, b_2dot)

num = -a'*b*(a_dot'*b + a'*b_dot);
num_dot = -a_dot'*b*(a_dot'*b + a'*b_dot) ...
    - a'*b_dot*(a_dot'*b + a'*b_dot) ...
    - a'*b*(a_2dot'*b + 2*a_dot'*b_dot + a'*b_2dot);

den = s(a, b);
den_dot = s_dot(a, b, a_dot, b_dot);

x = diff_num_den(num, num_dot, den, den_dot);

end