clear;
c=35e+7;
f=[20:1:150000];

warburg=1e11;
Rs=160000;
Re=10e8;
beta=0.87;
zcpe=c*(i*2*pi*f).^(-beta);
a=zcpe+Re;

a1=zcpe*Re;
z1=a1./a;
z=z1+Rs;
%chip_right_el26=zeros(3,17);
chip1_el2(1,:)=[130 270 520 1020 2020 4020 7000 10000 20000 40000 80000 150000];
chip1_el2(2,:)=[2320 1535 1130 609 414 239 161 127 73.3 43 24.2 15];
chip1_el2(3,:)=[-32.7 -31.1 -59 -67.5 -68.5 -677.8 -67.3 -68.8 -67.6 -68.7 -70.8 -66.7];
chip_right_el26(1,:)=[20 30 40 60 90 130 270 520 1020 2020 4020 7000 10000 20000 40000 80000 150000];
chip_right_el26(2,:)=[5550 3960 2970 1960 1506 1100 617 397 289 241 213 194 180 155 136 126 119];
chip_right_el26(3,:)=-[76 75 75 74 73 70 61 50 37 26 21 20 19 17 14 11 10];
fp=chip_right_el26(1,:);
modulp=chip_right_el26(2,:)*1e3;
fazap=chip_right_el26(3,:);

%imp=el_imp_mod1(c,beta,Re,warburg,Rs,f);

Ce=36e+7;
Rs=135000;
Re=10e8;
beta=0.86;

figure(1)
%subplot(1,2,1);
imp0=el_imp_mod3(Ce,beta,Re,Rs,0.5e-15,f);
%subplot(2,2,1);

fs=22;
a=loglog(f,abs(imp0)/1e3,fp,modulp/1e3,'b*');
axis([1e1 2e5 1e2 1e4])
set(a(1),'LineWidth',2);
set(a(1),'Color',[51/256,153/256,102/256]);
set(a(2),'MarkerSize',8)
grid on
h=gca;
set(h,'FontSize',fs);
set(h,'XTick',[1e1 1e2 1e3 1e4 1e5 1e6]);
set(h,'LineWidth',2);
%xlabel('f');
ylabel('|Z| [k\Omega]');

figure(2)
%subplot(1,2,2)
%subplot(2,2,2);
a=semilogx(f,angle(imp0)*180/pi,fp,fazap,'b*')
axis([1e1 2e5 -80 0])
set(a(1),'LineWidth',2);
set(a(1),'Color',[51/256,153/256,102/256]);
set(a(2),'MarkerSize',8)
grid on
h=gca;
set(h,'FontSize',fs);
get(h)
set(h,'XTick',[1e1 1e2 1e3 1e4 1e5 1e6]);
set(h,'LineWidth',2);
%xlabel('f');
ylabel('\phi [stopnie]');
%figure(3)

c=52e+7;
beta=0.95;
c2=9e+5;
beta2=0.15;

imp=el_imp_mod2(c,beta,Re,c2,beta2,f);
%subplot(2,2,3);
%loglog(f,abs(imp),fp,modulp);
grid on
%subplot(2,2,4);
%semilogx(f,angle(imp)*180/pi,fp,fazap);
grid on
model=struct('Ce',c,'beta',beta,'Re',Re,'Rs',Rs);
%m1=fit_el_mod1(chip1_el2',model);
%m1=fit_el_mod1(chip_right_el26',model);