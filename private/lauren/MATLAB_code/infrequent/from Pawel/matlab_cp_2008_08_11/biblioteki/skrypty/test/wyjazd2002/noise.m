run char1;
cd /home/pawel/pliki/nauka/neuro64/wyjazd2002/neuro64;

wzm1=wzm(2);
wzm2=wzm(4)

figure;

a=importdata('noise_14V.txt')*1500/1650;
b=importdata('noise_10V.txt')*1500/1650;
c=importdata('noise_08V.txt')*1500/1650;

t=[1:513];
plot(t,a(:,2),t,b(:,2),t,c(:,2));
axis([0 520 0 0.04]);
grid on;
legend('Vpolar=-1.4V','Vpolar=-1.0V','Vpolar=-0.8V');
xlabel('channel number');
ylabel('equivalent input noise level [mV]');
