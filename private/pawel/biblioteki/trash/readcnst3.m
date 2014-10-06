function s=readcnst3(name,header,nrchns,channels,samples);
a='bum';

lchannels=length(channels);

name0(1,1:7)=name(1,1:7);
name0(1,8:11)='conv';

name1(1,1:7)=name(1,1:7);
name1(1,8:11)='cnst';

fid0=fopen(name1,'r');
window=fread(fid0,1,'int32');
marg=fread(fid0,1,'int32');
%marg=1;
fclose(fid0);

l_dane=raportconv(name0,header,nrchns); %dlugosc kanalu w *conv

if samples(1)<1
    error('Zbyt male wartosci w "samples"');
end

if samples(2)>l_dane
    error('Zbyt duze wartosci w "samples"');
end

l_cnst=raportconv(name1,12,nrchns);
n_wezly=ceil(l_dane/window)+2;
cnst=zeros(1,n_wezly);  %"+1" - istotne ze wzgledu
%na przyjeta metode interpolacji - wartosci weltora x (argument
%funkcji "interp1") zmieniaja sie od -window/2 do length(sygnal)
%+window/2;

q=marg/window+0.5;
x_start=q-ceil(q);
x_stop=x_start+n_wezly-1;
%clear q;

x_wezly=[x_start:1:x_stop];
clear x_start;
clear x_stop;
x_wezly=floor(x_wezly.*window);

y_wezly=zeros(1,n_wezly);
%length(cnst)
y_start=ceil(q)+1;
y0=readconv(name1,12,nrchns,channels,[1 l_cnst]);
y_wezly(1,y_start:y_start+l_cnst-1)=y0;
y_wezly(1,1:(y_start-1))=y_wezly(1,y_start);
y_wezly(1,(y_start+l_cnst):(n_wezly))=y_wezly(1,y_start+l_cnst-1);
clear y0;

%y_wezly(
%indeks=find(x_wezly>=samples(1) & x_wezly<=samples(2))
%indeks1=indeks(1)-1;
%indeks2=max(indeks)+1;

%samples;

c=interp1(x_wezly,y_wezly,[samples(1):samples(2)],'nearest');

s=readconv(name0,header,nrchns,channels,samples)-c;   

a='bla';