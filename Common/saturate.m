function x_sat = saturate(x, x_min, x_max)

if x > x_max
    x_sat = x_max;
elseif x < x_min
    x_sat = x_min;
else
    x_sat = x;
end

end