n=64;
dt=0.001;

t=[-n*dt:dt:(n-1)*dt];

srednia=0;
sig=sqrt(n)/2*dt;

x=exp(-(t-srednia).^2/(2*sig^2));
subplot(2,1,1);
plot(x);
f=fft(x);
f=fftshift(f);
subplot(2,1,2);
plot(abs(f))

std(x)
std(abs(f))