cd /home2/pawel/oko/2000-12-12;

nrchns=65;
s1=239700;
dl=4000000;
samples=[s1 (s1+dl-1)];
channelsy=[7 17 28 29 34 65];
fs=20000;

threshold=300;

%channels=31;
name='Data009';

for kanal=1:6
	channels=channelsy(1,kanal);
	read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
	s=readcnst(read_param);
	time=[1:dl];

	detect_param=struct('prog',threshold,'histereza',50,'znak',-1);
	wynik=detect_f(read_param,detect_param);
	size(wynik)

	[odl,wys]=ampl_vs_opozn(read_param,detect_param);

	marg_left=20;
	marg_right=69;

	figure(4);
	%hold off;
	subplot(3,2,kanal);
	%clf;
	spikes=zeros(marg_left+marg_right+1);

	nadpr=4;
	prog=threshold;

	for i=1:length(wynik);
		start=wynik(1,i)-marg_left;
		stop=wynik(1,i)+marg_right;
		if (start>marg_left & stop<length(s)-marg_right)
			spike=s(start:stop);
			spikes(1,:)=spike;
			[spike0,filtr]=oversampling(spike,nadpr,11,0.9);
			szczyt=find(abs(spike0)==max(abs(spike0)));
			a=find(abs(spike0)>max(abs(spike0))/2);
			punkt0=a(1,1);
			szczyt=punkt0;
			if (szczyt>50 & szczyt<length(spike0)-270)
				spike1=spike0(1,szczyt-50:szczyt+269);
			end
			widmo1=abs(fft(spike1))/length(spike1);
			lw1=length(widmo1);
			f=[0:lw1-1]/lw1*fs*nadpr;
			%plot(spike);
			plot(f,widmo1);
			%loglog(f,widmo1);
			hold on;
		end
	end
	grid on;
	axis([0 round(fs/2) 0.01 200]);
	hold off;
end
