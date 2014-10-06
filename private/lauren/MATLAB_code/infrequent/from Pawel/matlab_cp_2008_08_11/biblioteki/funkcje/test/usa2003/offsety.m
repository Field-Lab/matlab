function y=offsety(filename,dlugosc,stim_channel);

start=1;

s=readconv(filename,206,65,stim_channel+1,[(start) (start+dlugosc-1)]);
sc=s-mean(s);
smin=max(max(abs(sc)));

detect=find(abs(sc)>smin*0.8);
poczatek=detect(1,1);

margines=poczatek+500;

for i=1:64;
	s=readconv(filename,206,65,i+1,[margines margines+100]);
	offsety(1,i)=mean(s);
end

y=offsety;
