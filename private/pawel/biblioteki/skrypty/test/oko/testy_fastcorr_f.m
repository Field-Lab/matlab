cd /home2/pawel/oko/2000-12-12/;

start=1;
dlugosc=100000;
%blok=2000000;
zasieg=[4000 4000];
wykladnik=15;

autokorelacje=zeros(64,8001);
korelacje65=autokorelacje;
korelacje27=autokorelacje;
korelacje28=autokorelacje;

clock
for i=1:64
    i
    read_param=struct('name','Data009','header',206,'nrchns',65,'channels',[i+1 i+1],'samples',[start start+dlugosc-1]);
    autokorelacje(i,:)=fastcorr_f(read_param,zasieg,wykladnik);
    
    read_param=struct('name','Data009','header',206,'nrchns',65,'channels',[65 i+1],'samples',[start start+dlugosc-1]);
    korelacje65(i,:)=fastcorr_f(read_param,zasieg,wykladnik);
    
    read_param=struct('name','Data009','header',206,'nrchns',65,'channels',[27 i+1],'samples',[start start+dlugosc-1]);
    korelacje27(i,:)=fastcorr_f(read_param,zasieg,wykladnik);
    
    read_param=struct('name','Data009','header',206,'nrchns',65,'channels',[28 i+1],'samples',[start start+dlugosc-1]);
    korelacje28(i,:)=fastcorr_f(read_param,zasieg,wykladnik);
end
clock

%save wyniki;
%exit;

%read_param=struct('name','Data009','header',206,'nrchns',65,'channels',[65 28],'samples',[start start+dlugosc-1]);

%clock
%z=fastcorr_f(read_param,zasieg,wykladnik);
%clock

%plot(z);