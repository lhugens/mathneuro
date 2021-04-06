function M = m(V)
    alpha_m = 0.055 .* (-27.01 - V) ./ (exp((-27.01 - V)/3.8) - 1);
    beta_m  = 0.94 .* exp((-63.01 - V)/17);
    M = alpha_m ./ (alpha_m + beta_m);
end
