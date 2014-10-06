%% Example found online
clear
a = 2;
dt = 0.05;
T = 100;
N = T/dt;
t = dt*(0:N-1)-50;
y = 2*(a/2/pi)./(t.^2+(a/2/pi)^2); % see website equation


df = 1/T;
f = df*(0:N-1);
Q = ceil((N+1)/2);
fQ = f(Q);
fc = f-fQ; % centering frequency

Y = (1/sqrt(2*pi)) * dt * fft(y); % see website equation
Y0 = sqrt(2*pi)*exp(-a*abs(fc)); % see website equation

plot(fc, Y0);
hold on
plot(fc, fftshift(abs(Y)),'r.'); % Slight discrepancy
axis([-10 10 0 3]);


%% 8 mm pupil, Artal & Navarro 1994
% I checked these again some actual PSFs from an older paper (1987?) and it
% looks right.

clear
df = 0.01; dt = 1 / df;
F = 200;
f = -F:df:F;

A = 0.53; B = 0.08; C = 0.11;
MTF = (1-C)*exp(-A*abs(f)) + C*exp(-B*abs(f));

x = dt/2*linspace(-1,1,length(f));
PSF = (1-C)*2*A./(A*A + 4*pi*pi*x.*x) + C*2*B./(B*B + 4*pi*pi*x.*x);

plot(x, df*abs(fftshift(fft(MTF))));
hold on;
plot(x, PSF, 'r.');

center = length(f) / 2;
interp1(PSF(center:end), x(center:end), max(PSF)/2) * 2


%% 6 mm pupil, Artal & Navarro 1994
clear
df = 0.01; dt = 1 / df;
F = 200;
f = -F:df:F;

A = 0.31; B = 0.06; C = 0.2;
MTF = (1-C)*exp(-A*abs(f)) + C*exp(-B*abs(f));

x = dt/2*linspace(-1,1,length(f));
PSF = (1-C)*2*A./(A*A + 4*pi*pi*x.*x) + C*2*B./(B*B + 4*pi*pi*x.*x);

plot(x, df*abs(fftshift(fft(MTF))));
hold on;
plot(x, PSF, 'r.');

center = length(f) / 2;
interp1(PSF(center:end), x(center:end), max(PSF)/2) * 2


%% 6 mm pupil, Guirao et al. 1999
clear
df = 0.01; dt = 1 / df;
F = 200;
f = -F:df:F;

a = 5.81; b = 21.68;
A = 1/a; B = 1/b; C = 1/4;
MTF = (1-C)*exp(-A*abs(f)) + C*exp(-B*abs(f));

x = dt/2*linspace(-1,1,length(f));
PSF = (1-C)*2*A./(A*A + 4*pi*pi*x.*x) + C*2*B./(B*B + 4*pi*pi*x.*x);

plot(x, df*abs(fftshift(fft(MTF))));
hold on;
plot(x, PSF, 'r.');

center = length(f) / 2;
interp1(PSF(center:end), x(center:end), max(PSF)/2) * 2