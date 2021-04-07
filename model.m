function dudt = model(t, u, I_app, P_max)
    % state variables
    V = u(1);
    C = u(2);

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
    %P_max = 0.002;

    % middle functions
    Csi     = (z .* F .* V) ./ (R .* T);
    I_drive = P_max .* z .* F .* Csi .* ((C - C_out .* exp(-Csi)) ./ (1 - exp(-Csi)));
    I_Ca_L  = m(V) .* h(C) .* I_drive;

    % time derivatives
    dudt    = zeros(2,1);
    dudt(1) = I_app(t) - g_L .* (V - E_L) - I_Ca_L;
    dudt(2) = -beta .* I_Ca_L - ((C - C_infty) ./ tau_C); 
end

% this is system input in MatCont
%{
z=2
F=96520
R=8313.4
T=273.15+25
C_out=2
g_L=0.05
E_L=-70
beyta=0.01
C_infty=0.0001
tau_C=200
P_max=0.002
ki=0.001;
h=ki/(ki+C)
alpha_m=0.055*(-27.01-V)/(exp((-27.01-V)/3.8)-1)
beta_m=0.94*exp((-63.01-V)/17)
m=alpha_m/(alpha_m+beta_m)
Csi=(z*F*V)/(R*T)
I_drive=P_max*z*F*Csi*((C-C_out*exp(-Csi))/(1-exp(-Csi)))
I_Ca_L=m*h*I_drive
V'=Iapp-g_L*(V-E_L)-I_Ca_L
C'=-beyta*I_Ca_L-((C-C_infty)/tau_C)
%}
