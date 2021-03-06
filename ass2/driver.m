clear all, close all, clc;

% exercise to run
% note: treat S2 of exercise 6 as exercise 7
VIEW_EXERCISE = 6;

%================================================================
if (VIEW_EXERCISE == 2)
end
%================================================================
if (VIEW_EXERCISE == 4)

    % m(V) - activation gating variable
    Vs = linspace(-100, 100, 1000);
    ms = m(Vs);

    % h(C) - inactivation gating variable
    Cs = linspace(-5, 5, 1000);
    hs = h(Cs);

    % m(V) plot and desired V value
    vplot = figure();
    subplot(2, 1, 1);
    plot(Vs, ms); 
    yline(0.9, 'r');
    xlabel('V'); ylabel('m'); grid on;
    disp(sprintf('V value for which m = 0.9: %.6e', fsolve(@(V) m(V) - 0.9, -20)));
    
    % h(C) plot
    subplot(2, 1, 2);
    plot(Cs, hs); 
    yline(0.1, 'r');
    xlabel('C'); ylabel('h'); grid on;
    disp(sprintf('C value for which h = 0.1: %.6e', fsolve(@(C) h(C) - 0.1, 0.001)));

end
%================================================================
if (VIEW_EXERCISE == 5)

    vplot = figure();
    plot_nullclines(1, true)
    xlabel('V'); ylabel('C'); grid on;
    
end
%================================================================
if (VIEW_EXERCISE == 6)

    % first, let's find the rest state for 
    % the case I_app = 0:

    % conditions
    P_max = 0.002;
    I_app = @(t) 0;

    % initial point
    u0(1) = 1.0;
    u0(2) = 1.0;

    % time evolution
    ts = [0 10000];
    dudt = @(t, u) model(t, u, I_app, P_max);
    [t, U] = ode45(dudt, ts, u0);
    Vs = U(:,1);
    Cs = U(:,2); 

    %disp(sprintf('Last coordinates of comet: V = %.6e; C = %.6e', Vs(end), Cs(end)));

    % approximate rest state for I_app = 0
    %u0(1) = -69.87100;
    %u0(2) =   0.01292298;

    % now, let's switch I_app to I at t = a

    CASE = 3;

    if     CASE == 1
        I = 4;
    elseif CASE == 2
        I = 1;
    elseif CASE == 3
        I = 1.2;
    end

    % conditions
    a = 3000;
    P_max = 0.002;
    I_app = @(t) I * heaviside(t - a);

    u0(1) = 1.0;
    u0(2) = 1.0;

    % time evolution
    ts = [0 10000];
    dudt = @(t, u) model(t, u, I_app, P_max);
    [t, U] = ode45(dudt, ts, u0);
    Vs = U(:,1);
    Cs = U(:,2); 

    % plot V(t)
    vplot = figure();
    subplot(2, 1, 1);
    plot(t, Vs);
    xlabel('t'); ylabel('V'); grid on;
    
    % plot C(t)
    subplot(2, 1, 2);
    plot(t, Cs); 
    xlabel('t'); ylabel('C'); grid on;

    % plot comet
    vplot = figure();
    plot_nullclines(I, false);
    xlabel('V'); ylabel('C'); grid on; hold on;
    xlim([-80 20]);
    ylim([0 2]);
    comet(Vs, Cs);

end

%================================================================
function plot_nullclines(I_app, signs)
    fp_V_dot = fimplicit(@(V,C) V_dot(V,C, I_app), [-100 100 -5 5]);
    fp_V_dot.Color = 'red';
    hold on;
    fp_V_dot = fimplicit(@(V,C) V_dot(V,C, 0), [-100 100 -5 5]);
    fp_V_dot.Color = 'black';
    hold on;
    fp_C_dot = fimplicit(@(V,C) C_dot(V,C), [-100 100 -5 5]);
    fp_C_dot.Color = 'blue';
    if (signs == true)
        Vs = linspace(-100, 100, 20);
        Cs = linspace(-5, 5, 20);
        r_v = 2.5;
        r_c = 0.25;
        for v = Vs
            for c = Cs
                v_dot = V_dot(v, c, I_app);
                c_dot = C_dot(v, c);
                if     (v_dot < 0) && (c_dot < 0)
                    theta = -3*pi/4;
                elseif (v_dot < 0) && (c_dot > 0)
                    theta =  3*pi/4;
                elseif (v_dot > 0) && (c_dot < 0)
                    theta =   -pi/4;
                elseif (v_dot > 0) && (c_dot > 0)
                    theta =    pi/4;
                end
                q = quiver(v, c, -r_v * cos(theta), -r_c * sin(theta), 'black', 'LineWidth', 0.8);
                q.ShowArrowHead = 'off';
                q.Annotation.LegendInformation.IconDisplayStyle = 'off';
                q.Marker = '.';
            end
        end
    end
end
