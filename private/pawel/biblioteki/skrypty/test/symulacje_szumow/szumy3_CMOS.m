%Szumy wzmacniacza CMOS s¹ modelowane dla tranzystorów PMOS, proces AMS
%0.35 um

%sta³e fizyczne:
kT=1.38e-3*300; %o dwadziescia rzedow wielkosci za duzo!!!

%Ogólne sta³e dla AMS 0.35 um:
Cox=4.5e-3; %unit: farad/square_meter
Muon=370e-4; %square meter/(Volt*sec) 
Muop=126e-4;
W=40e-6;
L=2e-6;
Id=3e-6;
%K=1e-25;

%Flicker noise constants:
KF=1.191e-26;
AF=1.461;

%Szumy CMOS
%a) szum termiczy pary ró¿nicowej
gmn=(2*Muon*Cox)^(1/2)*(W/L*Id)^(1/2);
gmp=(2*Muop*Cox)^(1/2)*(W/L*Id)^(1/2);
WhiteNoiseDensity=16/3*kT*(1e-20)/gmp
WhiteNoiseSigma=(WhiteNoiseDensity*1e4)^(1/2);

Flicker2=K/(Cox*W*L);

Rs=65e3;
%Rs=0;
Re=6e9;
Ce=1e-9;
%Ce=Ce*1000;

Cp=20e-12;
Rin=160e9;
omega_gr=1/(Cp*Rin);
f_gr=omega_gr/6.28
kT=1.38e-3*300; %o dwadziescia rzedow wielkosci za duzo!!!
 
f_step=0.1;
f=[10:f_step:5000];
omega=2*pi*f;
Ugm=ones(1,length(f))*WhiteNoiseDensity*1e20;
UFlicker=Flicker./f*1e20;

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
calk=Urin+Urs+Ure+Ugm+UFlicker;
loglog(f,Urin,f,Urs,f,Ure,f,Ugm,f,UFlicker);
axis([min(f) max(f) 1 1e9])
grid on;
legend('Urin','Urs','Ure','Ugm','Ufl');
nap_Rin=sqrt(sum(Urin))/1e10*sqrt(f_step)  %/1e10 - bo przedtem zbyt duza wartosc kT!!
nap_Rs=sqrt(sum(Urs))/1e10*sqrt(f_step) 
%nap_calk=sqrt(sum(calk))/1e10*sqrt(f_step)
nap_gm=sqrt(sum(Ugm))/1e10*sqrt(f_step)
nap_flick=sqrt(sum(UFlicker))/1e10*sqrt(f_step)
nap_calk=sqrt(sum(calk))/1e10*sqrt(f_step)
%sum(Ire)