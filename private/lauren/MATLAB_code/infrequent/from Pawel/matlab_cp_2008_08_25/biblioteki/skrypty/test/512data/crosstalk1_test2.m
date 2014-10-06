clear;

fs=20000;
wykladnik=14;

tic

cd I:\crosstalk_Sasha\data001;
filename='data001000.bin';
channel0=128;
f0=25;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data001.bin','wb')
fwrite(fid,data,'double');
fclose(fid);


cd I:\crosstalk_Sasha\data002;
filename='data002000.bin';
%channel0=128;
f0=50;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data002.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data003;
filename='data003000.bin';
%channel0=385;
f0=100;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data003.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data004;
filename='data004000.bin';
%channel0=385;
f0=200;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data004.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data005;
filename='data005000.bin';
%channel0=385;
f0=400;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data005.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data006;
filename='data006000.bin';
%channel0=385;
f0=800;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data006.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data007;
filename='data007000.bin';
%channel0=385;
f0=1600;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data007.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data008;
filename='data008000.bin';
%channel0=385;
f0=3200;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data008.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data009;
filename='data009000.bin';
%channel0=385;
f0=6400;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data009.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data010;
filename='data010000.bin';
%channel0=385;
f=9500;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data010.bin','wb')
fwrite(fid,data,'double');
fclose(fid);





cd I:\crosstalk_Sasha\data011;
filename='data011000.bin';
channel0=129;
f0=25;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data011.bin','wb')
fwrite(fid,data,'double');
fclose(fid);


cd I:\crosstalk_Sasha\data012;
filename='data012000.bin';
%channel0=128;
f0=50;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data012.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data013;
filename='data013000.bin';
%channel0=385;
f0=100;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data013.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data014;
filename='data014000.bin';
%channel0=385;
f0=200;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data014.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data015;
filename='data015000.bin';
%channel0=385;
f0=400;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data015.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data016;
filename='data016000.bin';
%channel0=385;
f0=800;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data016.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data017;
filename='data017000.bin';
%channel0=385;
f0=1600;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data017.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data018;
filename='data018000.bin';
%channel0=385;
f0=3200;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data018.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data019;
filename='data019000.bin';
%channel0=385;
f0=6400;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data019.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data020;
filename='data020000.bin';
%channel0=385;
f=9500;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data020.bin','wb')
fwrite(fid,data,'double');
fclose(fid);





cd I:\crosstalk_Sasha\data021;
filename='data021000.bin';
channel0=384;
f0=25;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data021.bin','wb')
fwrite(fid,data,'double');
fclose(fid);


cd I:\crosstalk_Sasha\data022;
filename='data022000.bin';
%channel0=128;
f0=50;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data022.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data023;
filename='data023000.bin';
%channel0=385;
f0=100;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data023.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data024;
filename='data024000.bin';
%channel0=385;
f0=200;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data024.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data025;
filename='data025000.bin';
%channel0=385;
f0=400;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data025.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data026;
filename='data026000.bin';
%channel0=385;
f0=800;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data026.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data027;
filename='data027000.bin';
%channel0=385;
f0=1600;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data027.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data028;
filename='data028000.bin';
%channel0=385;
f0=3200;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data028.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data029;
filename='data029000.bin';
%channel0=385;
f0=6400;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data029.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data030;
filename='data030000.bin';
%channel0=385;
f=9500;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data030.bin','wb')
fwrite(fid,data,'double');
fclose(fid);





cd I:\crosstalk_Sasha\data031;
filename='data031000.bin';
channel0=385;
f0=25;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data031.bin','wb')
fwrite(fid,data,'double');
fclose(fid);


cd I:\crosstalk_Sasha\data032;
filename='data032000.bin';
%channel0=128;
f0=50;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data032.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data033;
filename='data033000.bin';
%channel0=385;
f0=100;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data033.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data034;
filename='data034000.bin';
%channel0=385;
f0=200;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data034.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data035;
filename='data035000.bin';
%channel0=385;
f0=400;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data035.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data036;
filename='data036000.bin';
%channel0=385;
f0=800;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data036.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data037;
filename='data037000.bin';
%channel0=385;
f0=1600;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data037.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data038;
filename='data038000.bin';
%channel0=385;
f0=3200;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data038.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data039;
filename='data039000.bin';
%channel0=385;
f0=6400;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data039.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

cd I:\crosstalk_Sasha\data040;
filename='data040000.bin';
%channel0=385;
f=9500;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data040.bin','wb')
fwrite(fid,data,'double');
fclose(fid);





cd I:\crosstalk_Sasha\data041;
filename='data041000.bin';
channel0=512;
f0=25;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data041.bin','wb')
fwrite(fid,data,'double');
fclose(fid);


cd I:\crosstalk_Sasha\data042;
filename='data042000.bin';
%channel0=128;
f0=50;
[y,psdf,sigma]=crosstalk1(filename,fs,f0,channel0,wykladnik);
data=[psdf abs(y') sigma'];
cd H:/pliki/nauka/crosstalk/;
fid=fopen('data042.bin','wb')
fwrite(fid,data,'double');
fclose(fid);

toc
