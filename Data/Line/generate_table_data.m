clear all;
close all;

clc;

decoupled = load('decoupled.mat');
coupled = load('coupled.mat');
bres = load('brescianini.mat');
kooij = load('kooijman.mat');

[avg_ex, avg_eR, avg_f, converge_t] = extract_data(decoupled);
output_values("Proposed", avg_ex, avg_eR, avg_f, converge_t)

[avg_ex, avg_eR, avg_f, converge_t] = extract_data(coupled);
output_values("Coupled", avg_ex, avg_eR, avg_f, converge_t)

[avg_ex, avg_eR, avg_f, converge_t] = extract_data(kooij);
output_values("Kooij", avg_ex, avg_eR, avg_f, converge_t)

[avg_ex, avg_eR, avg_f, converge_t] = extract_data(bres);
output_values("Bres", avg_ex, avg_eR, avg_f, converge_t)

function output_values(name, avg_ex, avg_eR, avg_f, converge_t)
disp(name + " " + avg_ex + " & " + avg_eR + " & " + avg_f ...
    + " & " + converge_t);
end

function [avg_ex, avg_eR, avg_f, converge_t] = extract_data(data)

avg_ex = 0;
max_ex = 0;

avg_eR = 0;
max_eR = 0;

avg_f = 0;

converge_t = 0;
converge_ex = 0.02;
is_converged = false;

N = data.N;
for i = 1:N
    ex = data.e.x(:,i);
    eR = data.e.R(:,i);
    
    norm_ex = norm(ex, 2);
    norm_eR = norm(eR, 2);
    
    avg_ex = avg_ex + norm_ex;
    if norm_ex > max_ex
        max_ex = norm_ex;
    end
    
    avg_eR = avg_eR + norm_eR;
    if norm_eR > max_eR
        max_eR = norm_eR;
    end
    
    t = data.t(i);
    if norm(ex) < converge_ex
        if ~is_converged
            converge_t = t;
            is_converged = true;
        end
    end
    
    norm_f = norm(data.thr(:,i), 2);
    avg_f = avg_f + norm_f;
end

avg_ex = avg_ex / N;
avg_eR = avg_eR / N;
avg_f = avg_f / N;
end