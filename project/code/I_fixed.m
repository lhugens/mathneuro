function i = I_fixed(V, a)
    global C g_L E_L V_T delta_T tau_W;

    i = (a + g_L).*(V - E_L) - g_L .* delta_T .* exp((V - V_T) ./ delta_T);
end
