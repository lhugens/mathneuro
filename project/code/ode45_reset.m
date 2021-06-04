function [t U] = ode45_reset(dudt, ts, u0, reset)
    [t, U] = ode45(dudt, ts, u0);
end
