chns=[384 383 382];

frqs=[100 330 1000 3000 6000 9000];
V_05=zeros(length(chns),length(frqs));
V_08=zeros(length(chns),length(frqs));
V_14=zeros(length(chns),length(frqs));


V_08(1,:)=[1360 1680 1580 1075 500 255];
V_08(2,:)=[5 8.5 11 10 2 3];
V_08(3,:)=[10 7 2 6 2 2];

V_05(1,:)=[1390 1720 1640 970 340 250];
V_05(2,:)=[7.5 9 14 14 4 1.5];
V_05(3,:)=[19 15 7 8 4.5 3];

V_14(1,:)=[990 1630 1620 950 450 255];
V_14(2,:)=[9 6 16 12 2 3];
V_14(3,:)=[9 5 3 6 2 1.5];

figure(1);
subplot(1,2,1);
%clf;
semilogx(frqs,V_05(2,:)./V_05(1,:)*100,'k*-');
axis([100 10000 0 2]);
xlabel('frequency [Hz]');
ylabel('crosstalk [%]');
title('channel 63');
hold on;
semilogx(frqs,V_08(2,:)./V_08(1,:)*100,'k+-');
axis([100 10000 0 2]);
%hold on;
semilogx(frqs,V_14(2,:)./V_14(1,:)*100,'ko-');
axis([100 10000 0 2]);
%hold on;
legend('Vpolar=-0.5','Vpolar=-0.8','Vpolar=-1.4');
grid on;

%figure(2);
subplot(1,2,2);
%clf;
semilogx(frqs,V_05(3,:)./V_05(1,:)*100,'k*-');
axis([100 10000 0 2]);
xlabel('frequency [Hz]');
ylabel('crosstalk [%]');
title('channel 62');
hold on;
semilogx(frqs,V_08(3,:)./V_08(1,:)*100,'k+-');
axis([100 10000 0 2]);
%hold on;
semilogx(frqs,V_14(3,:)./V_14(1,:)*100,'ko-');
axis([100 10000 0 2]);
%hold on;
legend('Vpolar=-0.5','Vpolar=-0.8','Vpolar=-1.4');
grid on;