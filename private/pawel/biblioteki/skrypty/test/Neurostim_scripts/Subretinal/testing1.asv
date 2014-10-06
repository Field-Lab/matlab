cd C:\home\Pawel\nauka\StimFiles\subretinal\testy;
clear
file1='data.bin';
file2='data_d.bin';

h1=fopen(file1,'r','l');
clk1=fread(h1,'int32');
fclose(h1);

h1=fopen(file2,'r','l');
clk2=fread(h1,'int32');
fclose(h1);

max(abs(clk1-clk2))