cd H:\pliki\nauka\stymulacja\chip\testy\2006_05_19;
dac=[-127:127];
slope=[1:32];
offset=[1:32];

a=importdata('volt_r3_dcground_120k_channels_scan.dat');
figure(1);
for i=1:32
    s=a((i-1)*256+1:i*256,1);
    [slope1,lin_err1,dane,off]=dac_lin5(dac,s');
    slope(i)=slope1;
    offset(i)=off;
    subplot(4,8,i)
    plot(s)
    %hold on
    axis([0 256 -0.1 0.1]);
    grid on;
end

figure(11)
clear s;
for i=1:32
    s0=a((i-1)*256+1:i*256,1);    
    s=[s0(1:127,1)' (s0(128,1)+s0(129,1))/2 s0(130:256,1)'];
    a1=plot(dac,-s*1000,'k-')
    set(a1,'LineWidth',2);
    hold on
    axis([-140 140 -120 120]);
    grid on;
end
h=gca;
fs=18;
set(h,'FontSize',fs)
set(h,'LineWidth',2);
xlabel('DAC value');
ylabel('output voltage [mV]');
figure(2)

a1=plot(slope*1e3,'kd-');
xlabel('channel number');
ylabel('slope [mV/LSB]');
set(a1,'LineWidth',2);
set(a1,'MarkerSize',8)
axis([0 33 0.76 0.79]);
grid on;
%xlabel('channel number');
%ylabel('channel slope [mV/LSB]');
h=gca;
set(h,'FontSize',fs)
set(h,'LineWidth',2);
figure(3)
a2=plot(offset*1e3,'kd-');
set(a2,'LineWidth',2);
set(a2,'MarkerSize',8)
axis([0 33 -2 2]);
xlabel('channel number');
ylabel('offset [mV]');
h=gca;
set(h,'FontSize',fs)
set(h,'LineWidth',2);
grid on;