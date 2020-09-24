function desired = command_circle(t)
rad = 1;
w = 2*pi / 10;
height = 1;

desired.x = rad*[cos(w*t) - 1, sin(w*t), - height]';
desired.v = w*rad*[-sin(w*t), cos(w*t), 0]';
desired.x_2dot = w^2*rad*[-cos(w*t), -sin(w*t), 0]';
desired.x_3dot = w^3*rad*[sin(w*t), -cos(w*t), 0]';
desired.x_4dot = w^4*rad*[cos(w*t), sin(w*t), 0]';

w = 2*pi / 40;
desired.b1 = [cos(w*t), sin(w*t), 0]';
desired.b1_dot = w*[-sin(w*t), cos(w*t), 0]';
desired.b1_2dot = w^2*[-cos(w*t), -sin(w*t), 0]';
end