function W_dot = W_dot(V, W, a)
    global E_L tau_W;

    % state variables
    % V = u(1);
    % W = u(2);

    % time derivatives
    % dudt    = zeros(2,1);
    % dudt(1) = (I_app(t) - g_L * (V - E_L) + g_L * delta_T * exp((V - V_T)/delta_T) - W) / C;
    W_dot = (a * (V - E_L) - W) / tau_W;
end
