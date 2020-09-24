function thr = fM_to_thr(f, M, param)

d = param.d;
ctf = param.c_tf;

f_to_fM = [1, 1, 1, 1;
           0, -d, 0, d;
           d, 0, -d, 0;
           -ctf, ctf, -ctf, ctf];

thr = f_to_fM \ [f; M];
end