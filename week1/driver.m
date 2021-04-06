clear all, close all, clc;

%% WARM UP

p(1) =     1;           % Cm:   membrane capacitance [microFarads/cm^2]
p(2) =   120;           % gNa:  sodium conductance [milliSiemens/cm^3]
p(3) =    36;           % gK:   potassium conductance [milliSiemens/cm^3]
p(4) =   0.3;           % gL:   leak conductance [milliSiemens/cm^3]
p(5) =    50;           % eNa:  sodium Nernst potential [milliVolts]
p(6) =   -77;           % eK:   potassium Nernst potential [milliVolts]
p(7) = -54.4;           % eL:   leak reversal potential [milliVolts]
p(8) = 3^((20-6.3)/10); % phi:  temperature factor, see ET, equation 1.44.

u0(1) =   -60;  % Initial voltange v [milliVolts]
u0(2) =   0.0;  % Initial value for activation variable n
u0(3) =   0.0;  % Initial value for activation variable m
u0(4) =   0.0;  % Initial value for inactivation variable h

IAppFun = @(t) zeros(size(t));

dvdt = @(t, u) hh(t, u, p, IAppFun);

ts = [0 50];

[t, U] = ode45(dvdt, ts, u0);

subplot(3, 1, 1);
plot(t, U(:,1));
xlabel('t [ms]'); ylabel('v [mV]'); grid on;

subplot(3, 1, 2);
plot(t, U(:,2:4));
xlabel('t [ms]'); legend({'n(t)','m(t)','h(t)'}); grid on;

subplot(3, 1, 3);
plot(t, IAppFun(t));
xlabel('t [ms]'); ylabel('IApp [muA/cm^2]'); grid on;

%The system has a steady state whith values $(V, n, m, h) \approx (-65 mV, 0.3, 0.05, 0.6)
nrest = U(end, 2);
mrest = U(end, 3);
hrest = U(end, 4);

disp(sprintf('At rest, the proportion of open Na channels is %.6e', mrest^3*hrest))
disp(sprintf('At rest, the proportion of open K  channels is %.6e', nrest^4))


%% SMALL DEPOLARIZATION

IAppFun = @(t) (2 .* heaviside(t-16) .* heaviside(-t+18));

dvdt = @(t, u) hh(t, u, p, IAppFun);

ts = [0 50];

[t, U] = ode45(dvdt, ts, u0);

subplot(3, 1, 1);
plot(t, U(:,1));
xlabel('t [ms]'); ylabel('v [mV]'); grid on;

subplot(3, 1, 2);
plot(t, U(:,2:4));
xlabel('t [ms]'); legend({'n(t)','m(t)','h(t)'}); grid on;

subplot(3, 1, 3);
plot(t, IAppFun(t));
xlabel('t [ms]'); ylabel('IApp [muA/cm^2]'); grid on;

