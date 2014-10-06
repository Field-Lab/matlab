cd /mnt/win3/oko/2000-12-12/;

clear
start=1;
dlugosc=10000000;
blok=200000;
zasieg=[10000 10000];

read_param=struct('name','Data009','header',206,'nrchns',65,'channels',[65 28],'samples',[start start+dlugosc-1]);
detect_param=struct('prog',750,'histereza',50,'znak',-1);

tic
%s1=detect2_f(read_param,detect_param);
%s2=detect2_fb(read_param,detect_param,1000000);
%s=detect2_f(read_param,detect_param);
z=fastcorr_f(read_param,zasieg,16);
toc

tic
z=fastcorr_f(read_param,zasieg,17);
toc

tic
z=fastcorr_f(read_param,zasieg,18);
toc

tic
z=fastcorr_f(read_param,zasieg,19);
toc


%tic
%z1=fastcorr_fc(read_param,zasieg,18);
%toc
%size(s1)
%size(s2)

%r=readcnst5(read_param);