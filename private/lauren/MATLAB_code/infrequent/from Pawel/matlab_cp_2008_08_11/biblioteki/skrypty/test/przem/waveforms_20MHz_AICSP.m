cd H:\pliki\nauka\stymulacja\chip\testy\rozne;

a=importdata('bipolar_20MHz_16uA_19100_chip1_ch5.CSV');
b=importdata('bipolar_20MHz_64uA_4800_chip2_ch31.CSV');
%b=importdata('TEK00003.CSV');

%figure(2)
figure(1)
time=a(:,1)*1e6-3570.5;
h=plot(time,a(:,2)+0.3,'k-',time,b(:,2)-0.3,'k-')
set(h(1),'LineWidth',2)
set(h(2),'LineWidth',2)
%b=legend('19k\Omega','5k\Omega')
%set(b,'LineWidth',2)
axis([-3 33 -0.7 0.7])
h=gca
xlabel('time [\mus]');
ylabel('output voltage [V]');
set(h,'FontSize',18)
set(h,'LineWidth',2)

clear

a=importdata('bipolar_20MHz_64uAvolt_19100_chip1_ch5.CSV');
b=importdata('bipolar_20MHz_64uAvolt_4800_chip2_ch31.CSV');
%b=importdata('TEK00003.CSV');

%figure(2)
figure(2)
time=a(:,1)*1e6-23.15;
h=plot(time,a(:,2)+0.35,'k-',time,b(:,2)-0.35,'k-')
set(h(1),'LineWidth',2)
set(h(2),'LineWidth',2)
%b=legend('19k\Omega','5k\Omega')
%set(b,'LineWidth',2)
axis([-3 33 -0.8 0.8])
h=gca;
xlabel('time [\mus]');
ylabel('output voltage [V]');
set(h,'FontSize',18)
set(h,'LineWidth',2);