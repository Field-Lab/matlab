fidr=fopen('channel64.bin','r');
s1=fread(fidr,200000,'integer*2');
fclose(fidr);

dt=1/20000;
koniec=200000*dt;
t1=(0:dt:(koniec-dt))';
t2=1:200000;