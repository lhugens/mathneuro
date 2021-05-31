function Vd = V_discr(V)
    global C g_L E_L V_T delta_T tau_w a;
    Vd = (g_L ./ C).^2 .* (1 - exp((V - V_T) ./ delta_T)).^2 - (4.*a)./tau_w;
end
