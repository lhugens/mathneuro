function v_dot = V_dot(V, W, I_app, a)
    global C g_L E_L V_T delta_T tau_W;

    % state variables
    % V = u(1);
    % W = u(2);

    % time derivatives
    % dudt    = zeros(2,1);
    % dudt(1) = (I_app(t) - g_L * (V - E_L) + g_L * delta_T * exp((V - V_T)/delta_T) - W) / C;
    v_dot = (I_app - g_L * (V - E_L) + g_L * delta_T * exp((V - V_T)/delta_T) - W) / C;
    % dudt(2) = (a * (V - E_L) - W) / tau_W;
end
