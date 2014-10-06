frqs=[50 100 500 1000 2000];

set1=zeros(8,5);
set2=set1;
set3=set1;
set4=set1;

set1(1,:)=[650 1260 1660 1600 1340];
set1(2,:)=[520 1040 1510 1460 1210];
set1(3,:)=[630 1110 1400 1340 1100];
set1(4,:)=[630 1160 1500 1430 1180];
set1(5,:)=[510 1080 1500 1440 1190];
set1(6,:)=[520 1060 1440 1370 1130];
set1(7,:)=[585 1150 1510 1460 1200];
set1(8,:)=[540 1060 1420 1350 1120];

set2(1,:)=[720 1420 1850 1800 1550];
set2(2,:)=[530 1140 1700 1670 1420];
set2(3,:)=[670 1250 1600 1540 1300];
set2(4,:)=[690 1300 1680 1620 1380];
set2(5,:)=[550 1190 1660 1600 1370];
set2(6,:)=[575 1190 1640 1600 1350];
set2(7,:)=[630 1260 1670 1620 1390];
set2(8,:)=[565 1140 1550 1500 1270];

set3(1,:)=[420 1010 1650 1600 1340];
set3(2,:)=[340 840 1500 1470 1220];
set3(3,:)=[400 900 1400 1340 1100];
set3(4,:)=[420 940 1470 1420 1170];
set3(5,:)=[330 850 1450 1410 1180];
set3(6,:)=[350 870 1450 1400 1150];
set3(7,:)=[380 910 1480 1440 1180];
set3(8,:)=[350 830 1390 1340 1110];

set4(1,:)=[450 1110 1790 1770 1530];
set4(2,:)=[360 920 1650 1620 1380];
set4(3,:)=[440 980 1550 1510 1280];
set4(4,:)=[450 1040 1650 1610 1340];
set4(5,:)=[370 980 1700 1600 1410];
set4(6,:)=[360 940 1580 1540 1310];
set4(7,:)=[420 1010 1640 1610 1380];
set4(8,:)=[360 910 1530 1490 1270];

wzm=zeros(4,1)

a0=mean(set1);
wzm(1,1)=a0(3);
a0=mean(set2)
wzm(2,1)=a0(3);
a0=mean(set3)
wzm(3,1)=a0(3);
a0=mean(set4)
wzm(4,1)=a0(3);


%set1=set1/max(max(set1));
set1=20*log10(set1);

%set2=set2/max(max(set2));
set2=20*log10(set2);

%set3=set3/max(max(set3));
set3=20*log10(set3);

%set4=set4/max(max(set4));
set4=20*log10(set4);

figure(11);
subplot(2,2,1);
%clf;
for i=1:8
    semilogx(frqs,set1(i,:),'k*-');
    axis([10 10000 50 68]);
    title('settings1');
    xlabel('frequency [Hz]');
    ylabel('[gain [dB]');
    hold on;
end
grid on;
text(15,66,'Ihf=Ilf=24uA');
text(15,64,'Vpolar=-1.1V');

subplot(2,2,2);
%clf;
for i=1:8
    semilogx(frqs,set2(i,:),'k*-');
    axis([10 10000 50 68]);
    title('settings2');
    xlabel('frequency [Hz]');
    ylabel('[gain [dB]');
    hold on;
end
grid on;
text(15,66,'Ihf=Ilf=35uA');
text(15,64,'Vpolar=-1.1V');

subplot(2,2,3);
%clf;
for i=1:8
    semilogx(frqs,set3(i,:),'k*-');
    axis([10 10000 50 68]);
     title('settings3');
    xlabel('frequency [Hz]');
    ylabel('[gain [dB]');
    hold on;
end
grid on;
text(15,66,'Ihf=Ilf=24uA');
text(15,64,'Vpolar=-1.4V');

subplot(2,2,4);
%clf;
for i=1:8
    semilogx(frqs,set4(i,:),'k*-');
    axis([10 10000 50 68]);
    title('settings4');
    xlabel('frequency [Hz]');
    ylabel('[gain [dB]');
    hold on;
end
grid on;
text(15,66,'Ihf=Ilf=35uA');
text(15,64,'Vpolar=-1.4V');