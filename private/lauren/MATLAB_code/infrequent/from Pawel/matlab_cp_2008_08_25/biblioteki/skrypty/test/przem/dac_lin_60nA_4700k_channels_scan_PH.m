clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_19/
a=importdata('60nA_4700k_b4_c1_channels_scan.dat');

dac=[-127:127];
R=4684e3;;
fs=20;

figure(3)
clf;
slope=ones(1,32)*nan;
for i=[1:20 22:32]
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
figure(4)
plot(slope*1e9,'bd-');
grid on
%axis([0 33 32.8 34]);
h=gca;
set(h,'FontSize',fs);
a=xlabel('channel number');
set(a,'FontSize',fs);
a=ylabel('slope [nA/DAC step]');
set(a,'FontSize',fs);
rozrzut=std(slope([1:20 22:32]))/mean(slope([1:20 22:32]))*100