clear all, close all, clc;

global C g_L E_L V_T delta_T tau_w a;

C       = 281;
g_L     = 20;
E_L     = -70.6;
V_T     = -50.4;
delta_T = 2;
tau_w   = 144;
a       = 4;

% exercise to run
VIEW_PART = 3;

%================================================================
if (VIEW_PART == 1)
    plot_nullclines(0);
end

%================================================================
if (VIEW_PART == 2)
    I = 0;
    a = 3000;
    I_app = @(t) 0;

    u0(1) = -80.0;
    u0(2) = 0.0;

    % time evolution
    ts = [0 10000];
    dudt = @(t, u) model(t, u, I_app);
    [t, U] = ode45(dudt, ts, u0);
    Vs = U(:,1);
    ws = U(:,2); 

    % plot comet
    vplot = figure();
    plot_nullclines(I);
    xlabel('V'); ylabel('w'); grid on; hold on;
    comet(Vs, ws);
end

%================================================================
if (VIEW_PART == 3)
    Vs  = linspace(-100, -20, 1000);
    re1 = zeros(length(Vs));
    re2 = zeros(length(Vs));
    im1 = zeros(length(Vs));
    im2 = zeros(length(Vs));
    a = 1;
    c = a / tau_w;
    for i=1:length(Vs)
        b = (g_L/C) * (1 - exp((Vs(i) - V_T)/delta_T));
        p = [a b c];
        r = roots(p);
        re1(i) = real(r(1));
        re2(i) = real(r(2));
        im1(i) = imag(r(1));
        im2(i) = imag(r(2));
    end
    vplot = figure();
    subplot(2, 2, 1);
    plot(Vs, re1); 
    subplot(2, 2, 2);
    plot(Vs, re2); 
    subplot(2, 2, 3);
    plot(Vs, im1); 
    subplot(2, 2, 4);
    plot(Vs, im2); 

end

%================================================================
if (VIEW_PART == 4)
    Vs  = linspace(-100, 0, 1000);
    Vds = V_discr(Vs);
    
    vplot = figure();
    subplot(2, 1, 1);
    plot(Vs, Vds); 
    xlabel('V'); ylabel('V_discriminant'); grid on;
    V_discr_root_p = V_T + delta_T * log(1 + sqrt((4*a*C^2)/(tau_w*g_L^2)));
    V_discr_root_m = V_T + delta_T * log(1 - sqrt((4*a*C^2)/(tau_w*g_L^2)));
    
    disp('V_discr_root_p =');
    disp(V_discr_root_p);
    disp('V_discr_root_m =');
    disp(V_discr_root_m);
end


%================================================================
function plot_nullclines(I_app)
    fp_V_dot = fimplicit(@(V,w) V_dot(V, w, I_app), [-100 30 -1000 1000]);
    hold on;
    fp_w_dot = fimplicit(@(V,w) w_dot(V, w), [-100 30 -1000 1000]);
    hold on;
end
