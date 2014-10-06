cd /home2/pawel/usa_2004/data/March_19_Monkey;
filename='read20_stim19_at_threshold.xls';
%filename='no_cadmium_rec15_stim41.xls';
start=0;
length=40000;

y=readxls(filename,start,length);
plot(y);
grid on;

for i=1:0
	i
	start=(i-1)*length+1;
	y=readxls(filename,start,length);
	plot(y);
	grid on;	
	refresh(1);
	pause;
end


