function num_den = diff2_num_den(num, num_den, num_2dot, ...
    den, den_dot, den_2dot)

num = den^2*(den*num_2dot - num*den_2dot) ...
    - (den*num_den - num*den_dot)*2*den*den_dot;
den = den^4;

num_den = num / den;
end