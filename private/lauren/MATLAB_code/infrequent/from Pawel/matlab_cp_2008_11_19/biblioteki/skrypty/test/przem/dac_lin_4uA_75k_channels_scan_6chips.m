clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_13/
a=importdata('4uA_75k_b4_c1_channels_scan.dat');

cd H:/pliki/nauka/stymulacja/chip/testy/2006_07_12/
%a=importdata('4uA_75k_channels_scan_board3_chip_right.dat');
a=importdata('60nA_4700k_channels_scan_board1_chip_right.dat');

dac=[-127:127];
R=74900;
fs=20;

figure(5)
clf;
slope=ones(1,32)*nan;
for i=[1:32]
    s=a((i-1)*256+1:i*256)';
    [slope1,lin_err]=dac_lin5(dac,s);
    slope(i)=slope1/R;
    %figure(2)
    subplot(4,8,i);
    plot(dac,lin_err,'bd-');
    grid on;
    axis([-150 150 -2 2]);
    h=gca;
    set(h,'FontSize',14);
    set(h,'LineWidth',2)
    
    set(h,'YTick',[-2 -1 0 1 2]);
    if i~=25
        set(h,'XTickLabel',[]);
        set(h,'YTickLabel',[]);
    end
    if i==25
        set(h,'XTickLabel',[]);
        a1=xlabel('DAC');
        set(a1,'FontSize',16);
        set(h,'YTick',[-2 -1 0 1 2]);
        a1=ylabel('lin. error [LSB]');
        set(a1,'FontSize',16);
    end
end
figure(7)
a1=plot(slope*1e9,'bd-');
set(a1,'LineWidth',2);
set(a1,'MarkerSize',8)
grid on
%axis([0 35 32.8 35]);
h=gca;
set(h,'FontSize',fs);
set(h,'LineWidth',2)
a=xlabel('channel number');
set(a,'FontSize',fs);
a=ylabel('slope [nA/DAC step]');
set(a,'FontSize',fs);
