%%  Code from hw2
%%  Geoffrey Weiner, 10/14/2012
%%  Hodgkin-Huxley model
clear; close all
% Constants
C_m  =   1.0; % membrane capacitance, in uF/cm^2
g_Na = 120.0; % maximum conducances, in mS/cm^2
g_K  =  36.0;
g_L  =   0.3;
E_Na =  45.0; % Nernst reversal potentials, in mV
E_K  = -82.0;
E_L  = -59.387;

% Channel gating kinetics
% Functions of membrane voltage
alpha_m = @(V) 0.1.*(V+45.0)./(1.0 - exp(-(V+45.0) ./ 10.0));
beta_m  = @(V) 4.0.*exp(-(V+70.0) ./ 18.0);
alpha_h = @(V) 0.07.*exp(-(V+70.0) ./ 20.0);
beta_h  = @(V) 1.0./(1.0 + exp(-(V+40.0) ./ 10.0));
alpha_n = @(V) 0.01.*(V+60.0)./(1.0 - exp(-(V+60.0) ./ 10.0));
beta_n  = @(V) 0.125.*exp(-(V+70) ./ 80.0);

% Membrane currents (in uA/cm^2)
I_Na = @(V,m,h) g_Na .* m.^3 .* h .* (V - E_Na);
I_K  = @(V, n)  g_K  .* n.^4      .* (V - E_K);
I_L  = @(V)     g_L               .* (V - E_L);

%%  Problem 1 - plotting rate functions

Vsweep = -90 : 0.25 : 70;

figure();

plot(Vsweep, alpha_m(Vsweep), 'r-');
hold on;
plot(Vsweep, beta_m(Vsweep),  'r--');
plot(Vsweep, alpha_h(Vsweep), 'g-');
plot(Vsweep, beta_h(Vsweep),  'g--');
plot(Vsweep, alpha_n(Vsweep), 'b-');
plot(Vsweep, beta_n(Vsweep),  'b--');
hold off;

title('rate function plots');
xlabel('V (mV)');
ylabel('rate');
xlim([Vsweep(1) Vsweep(end)]);
legend('\alpha_m', '\beta_m', '\alpha_h', '\beta_h', '\alpha_n', '\beta_n', 'Location', 'SouthEastOutside');

%%  Problem 2 - leak current only, no channels

% injected current
% step up 10 uA/cm^2 every 100ms
I_inj = @(t) 10 * floor(t / 100)

dVmdt_leak = @(t, V) (I_inj(t) - I_L(V)) / C_m;
[t_leak, V_leak] = ode45(dVmdt_leak, [0 500], E_L);

figure();

subplot(2,1,1);
plot(t_leak, V_leak, 'k');
title('leaky current only, no channels');
ylabel('V (mV)');

subplot(2,1,2);
plot(t_leak, I_inj(t_leak), 'k');
xlabel('t (ms)');
ylabel('I_{inj} (\mu{A}/cm^2)');
ylim([-1 max(I_inj(t_leak))+1]);

%%  Problem 3 - full HH model

% injected current
% step up 10 uA/cm^2 every 100ms
I_inj = @(t) 10 * floor(t / 100);

% X = [V, m, h, n]
dVmdt = @(t, X) [
    (I_inj(t) - I_Na(X(1), X(2), X(3)) - I_K(X(1), X(4)) - I_L(X(1))) / C_m; ... % dV / dt
    alpha_m(X(1))*(1.0-X(2)) - beta_m(X(1))*X(2); ... % dm / dt
    alpha_h(X(1))*(1.0-X(3)) - beta_h(X(1))*X(3); ... % dh / dt
    alpha_n(X(1))*(1.0-X(4)) - beta_n(X(1))*X(4); ... % dn / dt
    ];

[t, X] = ode23(dVmdt, [0 500], [-70; 0.05; 0.6; 0.32]);

V = X(:,1);
m = X(:,2);
h = X(:,3);
n = X(:,4);

figure();
subplot(3,1,1);
plot(t, V, 'k');
title('full HH model')
ylabel('V (mV)');

subplot(3,1,2);
plot(t, m, 'r');
hold on;
plot(t, h, 'g');
plot(t, n, 'b');
hold off;
ylabel('Gating Value');
legend('m','h','n', 'Location', 'SouthWest');

subplot(3,1,3);
plot(t, I_inj(t), 'k');
xlabel('t (ms)');
ylabel('I_{inj} (\mu{A}/cm^2)');
ylim([-1 max(I_inj(t))+1]);

% injected current up to 200 uA/cm^2
% step up 50 uA/cm^2 every 100ms
I_inj = @(t) 50 * floor(t / 100)
dVmdt = @(t, X) [
    (I_inj(t) - I_Na(X(1), X(2), X(3)) - I_K(X(1), X(4)) - I_L(X(1))) / C_m; ... % dV / dt
    alpha_m(X(1))*(1.0-X(2)) - beta_m(X(1))*X(2); ... % dm / dt
    alpha_h(X(1))*(1.0-X(3)) - beta_h(X(1))*X(3); ... % dh / dt
    alpha_n(X(1))*(1.0-X(4)) - beta_n(X(1))*X(4); ... % dn / dt
    ];
[t, X] = ode23(dVmdt, [0 500], [-70; 0.05; 0.6; 0.32]);

V = X(:,1);
m = X(:,2);
h = X(:,3);
n = X(:,4);

figure();
subplot(3,1,1);
plot(t, V, 'k');
title('full HH model')
ylabel('V (mV)');

subplot(3,1,2);
plot(t, m, 'r');
hold on;
plot(t, h, 'g');
plot(t, n, 'b');
hold off;
ylabel('Gating Value');
legend('m','h','n', 'Location', 'SouthWest');

subplot(3,1,3);
plot(t, I_inj(t), 'k');
xlabel('t (ms)');
ylabel('I_{inj} (\mu{A}/cm^2)');
ylim([-1 max(I_inj(t))+1]);


%%  Problem 4 - threshold potential











