function [x, v, R, W, bar_theta_x, bar_theta_R] = split_to_states(X)

x = X(1:3);
v = X(4:6);
W = X(7:9);
R = reshape(X(10:18), 3, 3);
bar_theta_x = X(19:21);
bar_theta_R = X(22:24);

end