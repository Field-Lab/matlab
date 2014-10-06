function []=doubleaxis(dac,dane0,res)

R=res;

dane=[dane0(1,1:127) (dane0(1,128)+dane0(1,129))/2 dane0(1,130:256)];

mnoznik=1e6;

figure(1);
[a,b,c]=plotyy(dac,dane/R*mnoznik,dac,dane)
set(a(2),'Ylim',[-2.5 2.5]);
x=get(a(2),'Ylim')
set(a(1),'Ylim',x/R*mnoznik)
grid on

title('curr r1 19200k b4 chip1 data figure','FontSize',22,'FontWeight','demi');
xlabel('DAC (number)','FontSize',22,'FontWeight','demi');
grid on;

set(a,'FontSize',16);

set(get(a(1),'Ylabel'),'String','Current [uA]','FontSize',22);
set(get(a(2),'Ylabel'),'String','Voltage [V]','FontSize',22);
set(a(1),'YColor',[0,0,0]);
set(gca,'YTick',-0.15:0.05:0.15);
set(gca,'LineWidth',1.5);
set(a(2),'YColor',[0,0,0]);
%set(c,'FontSize',18)
get(b);
set(b,'Color',[0,0,1],'Marker','d');
set(c,'Color',[0,0,1],'Marker','d');
%set(b,'XMinorGrid','on')
