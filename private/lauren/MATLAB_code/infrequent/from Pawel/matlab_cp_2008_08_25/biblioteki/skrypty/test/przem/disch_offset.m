cd H:\pliki\nauka\stymulacja\chip\testy\2006_05_17;

figure(1)
for i=1:32
    filename=['probe_r4_120k_b' num2str(i-1) '_c1_dc400neg.dat'];
    a=importdata(filename);
    v0=a(2:2:82,1);
    v1=a(3:2:83);
    %subplot(4,8,i);
    plot(v0,v1,'bd-');
    hold on;
    %axis([1 41 -2e-3 2e-3]);
    grid on;
    m(i)=mean(v1-v0);
end
    
figure(2);
a1=plot(m*1e3,'bd-');
set(a1,'LineWidth',2);
set(a1,'MarkerSize',8)

a=xlabel('N');
set(a,'FontSize',18);
a=ylabel('offset [mV]');
set(a,'FontSize',18);
axis([0 33 -2.5 0.5])
grid on;
h=gca;
set(h,'FontSize',18);
set(h,'LineWidth',2);