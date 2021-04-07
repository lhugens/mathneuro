function H = h(C)
    ki = 0.001;
    H = ki ./ (ki + C);
end
