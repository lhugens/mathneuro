clear all, close all, clc;

global C g_L E_L V_T delta_T tau_W a;

g_L     = 20;
E_L     = -70.6;
V_T     = -50.4;
delta_T = 2;
tau_W   = 144;
C       = 281;

% type I
%a = 0.2*g_L;
%C = 3*tau_W*g_L;

% type II
%a = 3*g_L;
%C = 0.5*tau_W*g_L;

% exercise to run
VIEW_PART = 5;

% GENERAL VIEW OF NULLCLINE BEHAVIOUR
%================================================================
if (VIEW_PART == 1)
    a = 4;
    I = 400;
    vplot = figure();
    plot_nullclines(I, a);
    xlabel('V'); ylabel('W'); grid on;
    title(sprintf('a = %d ; I = %d', a, I))

    a = 100;
    I = 2400;
    vplot = figure();
    plot_nullclines(I, a);
    xlabel('V'); ylabel('W'); grid on;
    title(sprintf('a = %d ; I = %d', a, I))

end

% PLOT V_FIXED VS I_APP
%================================================================
if (VIEW_PART == 2)
    a = 4;
    vplot = figure();
    I_turn = plot_bifurcation(a);
    xlabel('I'); ylabel('V'); grid on;
    title(sprintf('a = %d ; I_{rh} = %0.2f', a, I_turn))

    a = 100;
    vplot = figure();
    I_turn = plot_bifurcation(a);
    xlabel('I'); ylabel('V'); grid on;
    title(sprintf('a = %d ; I_{rh} = %0.2f', a, I_turn))
end

% PLOT EIGENVALUES OF JACOBIAN
%================================================================
if (VIEW_PART == 3)
    Vs = linspace(-100, -40, 1000);
    [r1 r2 i1 i2] = J_eigen(Vs, 4);

    vplot = figure();
    subplot(2, 2, 1);
    plot(Vs, r1);
    subplot(2, 2, 2);
    plot(Vs, r2);
    subplot(2, 2, 3);
    plot(Vs, i1);
    subplot(2, 2, 4);
    plot(Vs, i2);

    idx = find(r2>0);
    V_turn = min(Vs(idx));
end

% SHOW COMET DURING BIFURCATION TYPE I - SADDLE-NODE
%================================================================
if (VIEW_PART == 4)
    a = 4;
    I = 401;

    tmax = 100000;
    Imin = 400;
    Imax = 500;
    m = (Imax - Imin)/tmax;
    b = Imin;
    I_app = @(t) m*t + b;

    u0(1) = -54.9;
    u0(2) = 62.7;

    % time evolution
    ts = [0 tmax];
    dudt = @(t, u) model(t, u, I_app, a);
    [t, U] = ode45(dudt, ts, u0);
    Vs = U(:,1);
    Ws = U(:,2); 

    % plot comet
    %vplot = figure();
    %plot_nullclines(I, a);
    %xlabel('V'); ylabel('W'); grid on; hold on;
    %comet(Vs, Ws);

    %vplot = figure();
    %subplot(2, 1, 1);
    %plot(t, Vs);
    %subplot(2, 1, 2);
    %plot(t, I_app(t));
end

% SHOW COMET DURING BIFURCATION TYPE II - HOPF
%================================================================
if (VIEW_PART == 5)
    a = 100;
    I = 2000;

    tmax = 1000000;
    Imin = 2000;
    Imax = 2600;
    m = (Imax - Imin)/tmax;
    b = Imin;
    I_app = @(t) m*t + b;

    u0(1) = -59.3;
    u0(2) = 2173.81;

    % time evolution
    ts = [0 tmax];
    dudt = @(t, u) model(t, u, I_app, a);
    [t, U] = ode45(dudt, ts, u0);
    Vs = U(:,1);
    Ws = U(:,2); 

    % plot comet
    vplot = figure();
    plot_nullclines(I, a);
    xlabel('V'); ylabel('W'); grid on; hold on;
    comet(Vs, Ws);

    vplot = figure();
    subplot(2, 1, 1);
    plot(t, Vs);
    subplot(2, 1, 2);
    plot(t, I_app(t));
end

%================================================================
function plot_nullclines(I_app, a)
    global E_L;
    Vmin = -100; 
    Vmax = 30;
    Wmin = a*(Vmin - E_L);
    Wmax = a*(Vmax - E_L);
    fp_V_dot = fimplicit(@(V,W) V_dot(V, W, I_app, a), [Vmin Vmax Wmin Wmax]);
    hold on;
    fp_W_dot = fimplicit(@(V,W) W_dot(V, W, a), [Vmin Vmax Wmin Wmax]);
end

function plot_V_fixed_vs_I(a)
    Vs = linspace(-100, 30, 1000);
    Is = I_fixed(Vs, a);
    plot(Is, Vs);
    hold on;

    Vs1 = linspace(-100, -40, 1000);
    [r1 r2 i1 i2] = J_eigen(Vs1, a);
    idx = find(r2>0);
    V_turn = min(Vs1(idx));
    idx = find(Vs < V_turn);
    V_unstable = Vs(idx);
    plot(I_fixed(V_unstable, a), V_unstable);

    xlim([-1000 max(Is)+100]);
end

function I_turn = plot_bifurcation(a)
    plot_V_fixed_vs_I(a)
    hold on;

    Vs = linspace(-100, -40, 1000);
    [r1 r2 i1 i2] = J_eigen(Vs, a);
    idx = find(r2>0);
    V_turn = min(Vs(idx));
    I_turn = I_fixed(V_turn, a);
    scatter(I_turn, V_turn, 'MarkerFaceColor', [.75 0 .75], 'MarkerEdgeColor', [.75 0 .75]); 
end
