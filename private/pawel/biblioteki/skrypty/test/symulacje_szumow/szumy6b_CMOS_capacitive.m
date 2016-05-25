clear
%Szumy wzmacniacza CMOS s¹ modelowane dla tranzystorów PMOS, proces AMS
%0.35 um. Tranzystory w parze ró¿nicowej pracuj¹w s³abej inwersji
%(2014-01-12).

%sta³e fizyczne:     
kT=1.38e-3*300; %o dwadziescia rzedow wielkosci za duzo!!! 300 to temepratura w Kelvinach
Ut=26e-3;

%Ogólne sta³e dla AMS 0.35 um:
Cox=4.5e-3; %unit: farad/square_meter
Muon=370e-4; %square meter/(Volt*sec) 
Muop=126e-4;
Kappa=0.7;

%Wymiary i polaryzacja tranzystorów:
W=50e-6/3;
L=18e-6/3;
Id=1e-6;

%Flicker noise constants:
% typical:
Tran=33; % 33 - tranzystory 3.3 V, 50 - tranzystory 5.0 V
Mode=2; % 1 - typical parameters, 2 - worst case
if Tran==33
    if Mode==1
        KF=1.191e-26;
        AF=1.461;
    else
        KF=1.827e-26;
        AF=1.405;
    end
else
    if Mode==1
        KF=3.625e-26;
        AF=1.638;
    else
        KF=4.852e-26;
        AF=1.564;
    end        
end

% 1. Szumy CMOS
% a) szum termiczy pary ró¿nicowej
gmp=Kappa*Id/Ut
WhiteNoiseDensity=16/3*kT*(1e-20)/gmp; %dla dwóch tranzystorów wg Majidzadeh2011, wzór 19
WhiteNoiseSigma=(WhiteNoiseDensity*1e4)^(1/2);

% b) szum 1/f pary ró¿nicowej
Flicker=KF/(Cox*W*L)*Id^AF/gmp^2*2; %the "*2" term comes from two transistors! Czy to prawda tak¿e dla s³abej inwrsji??

%c) parametry elektrody
Rs=60e3;
Re=6e12;
Ce=1e-9;

% d) parametry sprzê¿enia zwrotnego
Cin=10e-12;
Cf=0.5e-12;
Rf=0.32e12;

omega_gr=1/(Cf*Rf);
f_gr=omega_gr/6.28
 
f_step=0.05;
f=[0.1:f_step:10000];
omega=2*pi*f;
Ugm=ones(1,length(f))*WhiteNoiseDensity*1e20;
UFlicker=Flicker./f*1e20;

Zf=Rf./(1+j*omega*Rf*Cf);
Irf=4*kT*1e-20/Rf;
Urf=Irf*abs(Zf).^2*1e20*(Cf/Cin)^2;

Irs=4*kT*1e-20/Rs;
Urs=ones(1,length(f))*Irs*Rs^2*1e20;

figure(2);
calk=Urf+Urs+Ugm+UFlicker;
loglog(f,Urf,f,Urs,f,Ugm,f,UFlicker);
axis([min(f) max(f) 1 1e10])
grid on;
legend('Urf','Urs','Ugm','Ufl');
nap_Rf=sqrt(sum(Urf))/1e10*sqrt(f_step)  %/1e10 - bo przedtem zbyt duza wartosc kT!!
nap_Rs=sqrt(sum(Urs))/1e10*sqrt(f_step) 
nap_gm=sqrt(sum(Ugm))/1e10*sqrt(f_step)
nap_flick=sqrt(sum(UFlicker))/1e10*sqrt(f_step)
nap_calk=sqrt(sum(calk))/1e10*sqrt(f_step)