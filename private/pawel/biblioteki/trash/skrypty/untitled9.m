x0=-10;
x1=10;

dx=0.1;
x0=x0-dx/2;
x1=x1+dx/2;

x=[x0:dx:x1];
I=(sin(x)./x).^2;

figure(4);
plot(x,I);

figure(5);
semilogy(x,I);
