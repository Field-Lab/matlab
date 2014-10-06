cd /home/pawel/pliki/189e/pliki/dane/chip2_plat;

a=importdata('stim1_chip2.dat');
for i=1:32
	subplot(4,8,i);
	plot(a(i,:));
end
