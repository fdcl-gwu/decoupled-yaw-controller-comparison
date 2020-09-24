function x = s_dot(a, b, a_dot, b_dot)

x = -a'*b*(a_dot'*b + a'*b_dot) / s(a, b);

end