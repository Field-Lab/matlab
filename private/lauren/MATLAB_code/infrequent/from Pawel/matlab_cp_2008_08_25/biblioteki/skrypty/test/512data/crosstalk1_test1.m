clear;

fs=20000;
f0=200;
wykladnik=14;

cd I:\crosstalk_Sasha\data034;
filename='data034000.bin';
channel0=385;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data034.bin','wb')
fwrite(fid,data,'double');
fclose(fid);


f=[0:2^wykladnik-1]/2^wykladnik*fs;
df=fs/2^wykladnik;

figure(1);
semilogy(abs(y));
grid on;

figure(2);
plot(sqrt(sum(psdf')*df/2));
grid on;

figure(3);
x=psdf(120,:);
loglog(f,x);