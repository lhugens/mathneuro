function [r1 r2 i1 i2] = J_eigen(V, a)
    global C g_L E_L V_T delta_T tau_W;

    roots1 = zeros(length(V), 1);
    roots2 = zeros(length(V), 1);

    m2 = -(1/C);
    m3 = a/tau_W;
    m4 = -(1/tau_W);

    for i=1:length(V)
        m1 = -(g_L/C) * (1 - exp((V(i)-V_T)/delta_T));
        c_b = -(m1 + m4);
        c_c = m1*m4 - m2*m3;
        r = roots([1 c_b c_c]);
        roots1(i) = r(1);
        roots2(i) = r(2);
    end

    r1 = real(roots1);
    r2 = real(roots2);
    i1 = imag(roots1);
    i2 = imag(roots2);
end
