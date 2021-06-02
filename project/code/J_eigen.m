function r = J_eigen(V)
    global C g_L E_L V_T delta_T tau_W a;

    m1 = (g_L/C) * (1 - exp((V-V_T)/delta_T));
    m2 = -(1/C);
    m3 = a/tau_W;
    m4 = -(1/tau_W);

    c_b = -(m1 + m4);
    c_c = m1*m4 - m2*m3;
    r = roots([1 c_b c_c]);
end
