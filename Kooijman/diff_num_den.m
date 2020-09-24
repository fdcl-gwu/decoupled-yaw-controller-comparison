function x_dot = diff_num_den(num, num_dot, den, den_dot)
x_dot = (den*num_dot - num*den_dot) / den^2;
end
