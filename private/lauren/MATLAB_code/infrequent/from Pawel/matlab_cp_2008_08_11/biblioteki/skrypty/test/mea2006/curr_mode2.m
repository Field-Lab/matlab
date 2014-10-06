cd H:/pliki/nauka/mea2006/pomiary/electrode4;

a0=importdata('curr_range2_100kOhm.CSV');
a1=importdata('curr_range2_el4.CSV');
fs=18;
figure(1);

subplot(2,2,1);
a=plot(a0(:,1)*1000-0.025,a0(:,2)*10+0.125)
get(a)
set(a,'LineWidth',1.5);
axis([-0.1 0.5 -0.7 0.7]);
h=ylabel('output current [\muA]');
set(h,'FontSize',fs);
h=gca;
set(h,'Box','off');
set(h,'FontSize',fs);
set(h,'LineWidth',1.5);

subplot(2,2,3);
a=plot(a1(:,1)*1000-0.025,a1(:,2));
set(a,'LineWidth',1.5);
axis([-0.1 0.5 -0.22 0.14]);
a=xlabel('time [ms]');
set(a,'FontSize',fs)
a=ylabel('electrode potential [V]');
set(a,'FontSize',fs)
h=gca;
set(h,'Box','off');
set(h,'FontSize',fs);
set(h,'LineWidth',1.5);

h=gcf;
set(h,'Color',[1 1 1]);

b0=importdata('volt_range4_100kOhm.CSV');
b1=importdata('volt_range4_el4.CSV');
%figure(2);

subplot(2,2,2);
a=plot(b0(:,1)*1000-0.025,b0(:,2));
set(a,'LineWidth',1.5);
axis([-0.1 0.5 -0.3 0.3]);
a=ylabel('output voltage [V]');
set(a,'FontSize',fs)
h=gca;
set(h,'Box','off');
set(h,'FontSize',fs);
set(h,'LineWidth',1.5);

%figure(2)
subplot(2,2,4);
a=plot(b1(:,1)*1000-0.025,b1(:,2));
set(a,'LineWidth',1.5);
axis([-0.1 0.5 -0.3 0.3]);
a=xlabel('time [ms]');
set(a,'FontSize',fs)
a=ylabel('electrode potential [V]');
set(a,'FontSize',fs)
h=gca;
set(h,'LineWidth',1.5);
set(h,'Box','off');
set(h,'FontSize',fs);

h=gcf;
set(h,'Color',[1 1 1]);

%figure(3)
%plot(a1(:,1)*1000-0.025,a1(:,2)+0.0128,b1(:,1)*1000-0.025,b1(:,2)*(-0.6)-0.007);