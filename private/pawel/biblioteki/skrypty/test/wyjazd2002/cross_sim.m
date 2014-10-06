c1=150e-12;
c2=50e-12;
c3=1e-12;

R1=[1e7 13e7 1000e7];

f=[1:1:10000];

j=sqrt(-1);
z1=(j*6.28*f*c1).^(-1);
z2=(j*6.28*f*c2).^(-1);
z3=(j*6.28*f*c3).^(-1);

figure(1);
clf;

for i=1:3
R=R1(1,i);
x0=z2.*(z1+R);
x1=z3.*(z1+z2+R);
x2=z2.*(z1+R);
x(i,:)=x0./(x1+x2);
end

x=x*100;

semilogx(f,x(1,:),'k-',f,x(2,:),'k-.',f,x(3,:),'k--');
grid on;
xlabel('frequency [Hz]');
ylabel('crasstalk [%]');
legend('Rin=10 MOhm','Rin=130 MOhm','Rin=10 GOhm',4);