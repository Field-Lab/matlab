clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_19/
a=importdata('60nA_4700k_channels_scan.dat');
dac=[-127:127];
R=4700e3;

figure(3)
clf;

slope=ones(1,32)*nan;
smin=slope;
for i=[1:20 22:32]
    s=a((i-1)*256+1:i*256)';
    [slope1,lin_err]=dac_lin3(dac,s);
    slope(i)=slope1/R;
    %figure(2)
    subplot(4,8,i);
    plot(dac,lin_err,'bd-');
    grid on;
    axis([-150 150 -2 2]);
    %pause(0.1);
    smin(i)=s(139);
end
figure(4)
plot(slope);
grid on
h=gca;
set(h,'FontSize',18);
figure(5);
plot(smin/R);