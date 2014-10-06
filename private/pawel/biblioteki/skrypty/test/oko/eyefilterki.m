n=1500;

fp=20000;
f1=15;
f2=3000;
f3=50;
f0=1;

r=eyefilter(n,f1,f2,f3,fp);
r1=inveyefilter(n,f3,f1,f2,fp,f0);

figure(1);

o=conv(r,r1);

subplot(2,2,1);
plot(r1);
%axis([1 4*n -0.0001 0.0001]);
grid on;

subplot(2,2,2);
semilogy(abs(o));
axis([0 4*n 1e-8 1]);
grid on;

z=zeros(1,8192);
z1=z;

z(1,1:(2*n+1))=r;
z1(1,1:(2*n+1))=r1;

f=fft(z);
f1=fft(z1);

subplot(2,2,3);
loglog(abs(f1));
grid on;

subplot(2,2,4);
plot(abs(f.*f1));
%axis([1 8192 0.98 1.02]);
grid on;

%figure(3);
%loglog(abs(f));
%grid on;

