clear
cd H:/pliki/nauka/stymulacja/chip/testy/CHANNELS_SCAN/
%a=importdata('250uA_1190_c1_channels_scan.dat');
%a=importdata('4uA_75k_c1_channels_scan.dat');
a=importdata('60nA_4700k_c1_channels_scan.dat');
dac=[-127:127];
R=1190;

figure(3)
clf;
slope=ones(1,32)*nan;
for i=[1:20 22:32]
    s=a((i-1)*256+1:i*256)';
    [slope1,lin_err]=dac_lin3(dac,s);
    slope(i)=slope1/R;
    %figure(2)
    subplot(4,8,i);
    plot(dac,lin_err,'bd-');
    grid on;
    axis([-150 150 -2 2]);
    pause(0.1);
end
figure(4)
plot(slope*127);
grid on
h=gca;
set(h,'FontSize',18);
title('Ch. gain vs ch. number r3 75k','FontSize',22,'FontWeight','demi');
xlabel('Channel number','FontSize',22,'FontWeight','demi');
ylabel('Channel gain','FontSize',22,'FontWeight','demi');
grid on;
h=gca;
set(h,'FontSize',16);
set(h,'LineWidth',1.5);