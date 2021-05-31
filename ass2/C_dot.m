function c_dot = C_dot(V, C)
    % state variables
    %V = u(1);
    %C = u(2);

    % constants
    z       = 2;
    F       = 96520;
    R       = 8313.4;
    T       = 273.15 + 25;
    C_out   = 2;
    g_L     = 0.05;
    E_L     = -70;
    beta    = 0.01;
    C_infty = 0.0001;
    tau_C   = 200;
    %I_app = 1;
    P_max = 0.002;

    % middle functions
    Csi     = (z .* F .* V) ./ (R .* T);
    I_drive = P_max .* z .* F .* Csi .* ((C - C_out .* exp(-Csi)) ./ (1 - exp(-Csi)));
    I_Ca_L  = m(V) .* h(C) .* I_drive;

    % time derivatives
    %dudt    = zeros(2,1);
    %dudt(1) = I_app(t) - g_L .* (V - E_L) - I_Ca_L;
    c_dot   = -beta .* I_Ca_L - ((C - C_infty) ./ tau_C); 
end
