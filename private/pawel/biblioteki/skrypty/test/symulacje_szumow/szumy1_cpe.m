Rs=250e3;
%Rs=0;
Re=600e6;
Ce=1.5e-9;
%Ce=Ce*1000;

Cp=15e-12;
Rin=500e6;
kT=1.38e-3*300; %o dwadziescia rzedow wielkosci za duzo!!!

f=[1:0.1:20];
omega=2*pi*f;

Zc=Re./(1+j*omega*Re*Ce)+1./(j*omega*Cp);

a=abs(Rs+Rin+Zc);
Ct=4*kT*1./(a.^2);
clear a;

Irs=Ct*Rs;

a=abs(Rs+Zc);
Irin=Ct.*a.^2/Rin;
clear a;

a=abs(1+j*omega*Re*Ce);
Ire=Ct*Re./a.^2;

Urin=Irin*Rin^2;
Urs=Irs*Rin^2;
Ure=Ire*Rin^2;

figure(1);
calk=Urin+Urs+Ure;
loglog(f,Urin,f,Urs,f,Ure,f,calk);
axis([1 max(f) 1 1e9])
grid on;
legend('Irin','Irs','Ire','calk.');
nap_Rin=sqrt(sum(Urin))/1e10  %/1e10 - bo przedtem zbyt duza wartosc kT!!
nap_Rs=sqrt(sum(Urs))/1e10
nap_calk=sqrt(sum(calk))/1e10
%sum(Ire)