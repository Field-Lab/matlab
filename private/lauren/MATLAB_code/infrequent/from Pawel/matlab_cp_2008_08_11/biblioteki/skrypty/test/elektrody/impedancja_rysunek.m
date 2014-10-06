clear;
c=35e+7;
f=[0.2:0.2:500000];

warburg=1e11;
Rs=160000;
Re=10e8;
beta=0.85;
zcpe=c*(i*2*pi*f).^(-beta);
a=zcpe+Re;

a1=zcpe*Re;
z1=a1./a;
z=z1+Rs;
%chip_right_el26=zeros(3,17);
Ce=25e+7;
Rs=150000;
Re1=2.5e9;
Re2=5e6;
Re3=1e7;
beta=0.82;

imp1=el_imp_mod1(Ce,beta,Re1,warburg,Rs,f);
%imp2=el_imp_mod1(Ce,beta,Re2,warburg,Rs,f);
%imp3=el_imp_mod1(Ce,beta,Re3,warburg,Rs,f);

%modul=abs(imp);
faza=angle(imp1);
figure(1)
%a=loglog(f,abs(imp1),f,abs(imp2),f,abs(imp3))
%axis([1e-2 1e7 1e4 1e8])
%set(a,'LineWidth',2)
%set(a,'Color',[51/256,153/256,102/256])
%grid on
figure(2)
semilogx(f,faza*180/pi);
%grid on
%model=struct('Ce',c,'beta',beta,'Re',Re,'Rs',Rs);
