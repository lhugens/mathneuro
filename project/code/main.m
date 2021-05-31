clear all, close all, clc;

global C g_L E_L V_T delta_T tau_w a;

C       = 281;
g_L     = 20;
E_L     = -70.6;
V_T     = -50.4;
delta_T = 2;
tau_w   = 144;
a       = 4;

Vs  = linspace(-46.8, -47.8, 1000);
Vds = V_discr(Vs);

vplot = figure();
subplot(2, 1, 1);
plot(Vs, Vds); 
xlabel('V'); ylabel('V_discriminant'); grid on;
V_discr_root_p = V_T + delta_T * log(1 + sqrt((4*a*C^2)/(tau_w*g_L^2)));
V_discr_root_m = V_T + delta_T * log(1 - sqrt((4*a*C^2)/(tau_w*g_L^2)));
disp(sprintf('V_discr_root_p = %.6e', V_discr_root_p));
disp(sprintf('V_discr_root_m = %.6e', V_discr_root_m));
