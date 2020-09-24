function [f_sat, M_sat] = saturate_fM(f, M, param)
thr = fM_to_thr(f, M, param);

max_f = 8;
min_f = 0.1;

for i = 1:4
    if thr(i) > max_f
        thr(i) = max_f;
    elseif thr(i) < min_f
        thr(i) = min_f;
    end
end

[f_sat, M_sat] = thr_to_fM(thr, param);
end