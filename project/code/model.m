function dudt = model(t, u, I_app, a)
    global C g_L E_L V_T delta_T tau_W;

    % state variables
    V = u(1);
    W = u(2);

    % time derivatives
    dudt    = zeros(2,1);
    dudt(1) = (I_app(t) - g_L * (V - E_L) + g_L * delta_T * exp((V - V_T)/delta_T) - W) / C;
    dudt(2) = (a * (V - E_L) - W) / tau_W;
end

%{
MatCont Input:

g_L=20
E_L=-70.6
V_T=-50.4
delta_T=2
tau_W=144
C=281
a=100
V'=(Iapp-g_L*(V-E_L)+g_L*delta_T*exp((V-V_T)/delta_T)-W)/C
W'=(a*(V-E_L)-W)/tau_W
%}

