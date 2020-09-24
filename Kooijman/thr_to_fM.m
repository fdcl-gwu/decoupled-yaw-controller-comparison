function [f, M] = thr_to_fM(thr, param)
d = param.d;
ctf = param.c_tf;

f_to_fM = [1, 1, 1, 1;
           0, -d, 0, d;
           d, 0, -d, 0;
           -ctf, ctf, -ctf, ctf];

fM = f_to_fM * thr;

f = fM(1);
M = fM(2:4);
end