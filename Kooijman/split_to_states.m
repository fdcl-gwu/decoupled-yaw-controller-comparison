function [x, v, R, W] = split_to_states(X)

x = X(1:3);
v = X(4:6);
R = reshape(X(7:15), 3, 3);
W = X(16:18);

end