function u_sat = saturate_u(u, a, b)

u_sat = zeros(3, 1);

for i = 1:3
    ui = u(i);
    ai = a(i);
    bi = b(i);
    
    e = 0.01;
    e_upper = (bi - ai) / 2;
    if e > e_upper
        e = e_upper;
    end
    
    if ai + e < ui && ui < bi - e
        u_sat(i) = u(i);
    elseif ui <= ai - e
        u_sat(i) = ai;
    elseif bi + e <= ui
        u_sat(i) = bi;
    elseif ai - e < ui && ui <= ai + e 
        u_sat(i) = ui + 1 / (4*e) * (ui - (ai + e))^2;
    elseif bi - e <= ui && ui < bi + e 
        u_sat(i) = ui - 1 / (4*e) * (ui - (bi - e))^2;
    end
end

end