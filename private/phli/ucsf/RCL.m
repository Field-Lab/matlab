R  = 3.68e9;
Rl = 4.04e9;
L  = 2.21e9;
C  = 6e-12;
w = [0 10.^(-1:0.1:2)];

f = w ./ 2 ./ pi;

% My calc
jw = 1i * w;
G = 1 / R;
Zl = L * jw;
Zll = Rl + Zl;
All = 1 ./ Zll;
Ac = C * jw;
Ztot = 1 ./ (G + All + Ac);


% Demontis et al. 1999
common = Rl^2 + w.^2 .* L^2;
denom = G + Rl ./ common;
num = w.*C - w.*L ./ common;
modulus = 1 ./ sqrt(denom.^2 + num.^2);
fase = -atan(num ./ denom);


%% Compare
subplot(3, 1, 1:2);
semilogx(f, modulus);
hold on
semilogx(f, abs(Ztot), 'r.');

subplot(3, 1, 3);
semilogx(f, fase);
hold on
semilogx(f, phase(Ztot), 'r.');