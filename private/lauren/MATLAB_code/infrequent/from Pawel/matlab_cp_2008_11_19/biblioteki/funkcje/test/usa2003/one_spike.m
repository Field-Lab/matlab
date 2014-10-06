function y=one_spike(filename,dlugosc,okres,stim_channel,channels,zakres);
%dlugosc, zakres - w ilosciach probek

y=zeros(length(channels),zakres);
start=1;

s=readconv(filename,206,65,stim_channel,[(start) (start+dlugosc-1)]);
sc=s-mean(s);
mean(s);
smin=max(max(abs(sc)));

detect=find(abs(sc)>smin*0.8);
poczatek=detect(1,1);

margines=poczatek+okres-50;

s=readconv(filename,206,65,stim_channel,[margines margines+okres]);
sc=s-mean(s(1,1:20));
mean(s(1,1:10))
smin=max(max(abs(sc)))

detect=find(abs(sc)>100);
poczatek=detect(1,1)

margines2=poczatek+okres-100;


for i=1:length(channels);
	s=readconv(filename,206,65,channels(i)+1,[margines+margines2 margines+margines2+zakres-1]);
	y(i,:)=s;%-mean(s(1,1:20));
end


