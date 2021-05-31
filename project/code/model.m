function dudt = model(t, u, I_app)
    global C g_L E_L V_T delta_T tau_w a;

    % state variables
    V = u(1);
    w = u(2);

    % time derivatives
    dudt    = zeros(2,1);
    %dudt(1) = (I_app(t) - g_L * (V - E_L) + g_L * delta_T * exp((V - V_T)/delta_T) - w) / C;
    dudt(1) = (I_app - g_L * (V - E_L) + g_L * delta_T * exp((V - V_T)/delta_T) - w) / C;
    dudt(2) = (a * (V - E_L) - w) / tau_w;
end
