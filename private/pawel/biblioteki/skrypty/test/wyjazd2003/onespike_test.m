cd /home2/pawel/oko/santa_cruz_2003/2003_10_20;
filename='ch28_25uA_01ms_Vpol-11_waterconv';

dlugosc=20000;
stim_channel=28;
channels=[28 35 45];
zakres=8800;
okres=4000;

time=[1:zakres]/20;

filename='ch28_25uA_10ms_Vpol-14_waterconv';
y1=one_spike(filename,dlugosc,okres,stim_channel,channels,zakres);

filename='ch28_25uA_10ms_Vpol-11_waterconv';
y2=one_spike(filename,dlugosc,okres,stim_channel,channels,zakres);

filename='ch28_25uA_10ms_Vpol-08_waterconv';
y3=one_spike(filename,dlugosc,okres,stim_channel,channels,zakres);

style1='b- ';
style2='g-- ';
style3='r:- ';

figure(4);
clf;
for i=1:3
	subplot(1,3,i);
	a1=plot(time,y1(i,:),'b-',time,y2(i,:),'g-',time,y3(i,:),'r-');
	grid on;
	%set(a1(2),'LineWidth',2);
	%set(a1(3),'LineWidth',2);

	axis([0 time(zakres) -1200 1200]);
	legend('R_{in}=10M\Omega','R_{in}=20M\Omega','R_{in}=150M\Omega');
	xlabel('time [miliseconds]');
	if (i==1) 
		ylabel('signal level [ADC units]');
	end
	title(['channel' num2str(channels(i))]);
end
