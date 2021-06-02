clear all, close all, clc;

global C g_L E_L V_T delta_T tau_W a;

g_L     = 20;
E_L     = -70.6;
V_T     = -50.4;
delta_T = 2;
tau_W   = 144;
a       = 4;
C       = 281;

% type I
%a = 0.2*g_L;
%C = 3*tau_W*g_L;

% type II
%a = 3*g_L;
%C = 0.5*tau_W*g_L;

% exercise to run
VIEW_PART = 2;

% GENERAL VIEW OF NULLCLINE BEHAVIOUR
%================================================================
if (VIEW_PART == 1)
    vplot = figure();
    subplot(1, 2, 1);
    plot_nullclines(0);
    subplot(1, 2, 2);
    plot_nullclines(500);
end

% PLOT V_FIXED VS I_APP
%================================================================
if (VIEW_PART == 2)
    vplot = figure();
    subplot(2, 1, 2);
    a = 100;
    xlabel('I'); ylabel('V'); grid on;
    plot_bifurcation();
    subplot(2, 1, 1);
    a = 4;
    xlabel('I'); ylabel('V'); grid on;
    plot_bifurcation();
end

% PLOT EIGENVALUES OF JACOBIAN
%================================================================
if (VIEW_PART == 3)
    Vs = linspace(-100, -40, 1000);
    [r1 r2 i1 i2] = J_eigen(Vs);

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


%================================================================
if (VIEW_PART == 10)
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
    Ws = U(:,2); 

    % plot comet
    vplot = figure();
    plot_nullclines(I);
    xlabel('V'); ylabel('W'); grid on; hold on;
    comet(Vs, Ws);
end

%================================================================
function plot_nullclines(I_app)
    fp_V_dot = fimplicit(@(V,W) V_dot(V, W, I_app), [-100 30 -1000 1000]);
    hold on;
    fp_W_dot = fimplicit(@(V,W) W_dot(V, W), [-100 30 -1000 1000]);
end

function plot_V_fixed_vs_I()
    Vs = linspace(-100, 30, 1000);
    Is = I_fixed(Vs);
    plot(Is, Vs);
    hold on;

    Vs1 = linspace(-100, -40, 1000);
    [r1 r2 i1 i2] = J_eigen(Vs1);
    idx = find(r2>0);
    V_turn = min(Vs1(idx));
    idx = find(Vs < V_turn);
    V_unstable = Vs(idx);
    plot(I_fixed(V_unstable), V_unstable);

    xlim([-1000 max(Is)+100]);
end

function plot_bifurcation()
    plot_V_fixed_vs_I()
    hold on;

    Vs = linspace(-100, -40, 1000);
    [r1 r2 i1 i2] = J_eigen(Vs);
    idx = find(r2>0);
    V_turn = min(Vs(idx));
    I_turn = I_fixed(V_turn);
    scatter(I_turn, V_turn, 'MarkerFaceColor', [.75 0 .75], 'MarkerEdgeColor', [.75 0 .75]); 
end
