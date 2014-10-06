Vt=26e-3*312/310;
beta1=0.45;
beta2=0.5;
beta3=0.55;
J0=25e-9;
dV=1e-4;
Vmax=0.5;
V=[-Vmax:dV:Vmax];
z=1;

J1=J0*(exp((1-beta1)*z*V/Vt)-exp(-beta1*z*V/Vt));
J2=J0*(exp((1-beta2)*z*V/Vt)-exp(-beta2*z*V/Vt));
J3=J0*(exp((1-beta3)*z*V/Vt)-exp(-beta3*z*V/Vt));

Fs=23;

figure(1)
a=plot(V*1000,J2)%,J1,V*1000,J2,V*1000,J3);
set(a,'LineWidth',2)
set(a,'Color',[51/256,153/256,102/256])
%axis([-400 400 -470*J0 470*J0])
h=gca;
set(h,'FontSize',Fs);
h=xlabel('Vnp [mV]');
set(h,'FontSize',Fs);
h=ylabel('J/J_{0}');
set(h,'FontSize',Fs);
grid on;

Rt=Vt/J0/z;
linia=V/Rt;
figure(2)
a2=plot(V*1000,J2,V*1000,linia);
a=a2(1)
set(a,'LineWidth',2)
set(a,'Color',[51/256,153/256,102/256])
a=a2(2)
set(a,'LineStyle','--');
set(a,'LineWidth',2)
set(a,'Color',[100/256,100/256,153/256])

%axis([-100 100 -7 7])
h=gca;
set(h,'FontSize',Fs);
h=xlabel('Vnp [mV]');
set(h,'FontSize',Fs);
h=ylabel('J/J_{0}');
set(h,'FontSize',Fs);
grid on