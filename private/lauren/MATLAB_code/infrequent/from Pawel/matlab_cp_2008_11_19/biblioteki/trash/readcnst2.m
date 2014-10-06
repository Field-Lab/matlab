function s=readcnst2(name,header,nrchns,channels,samples);

lchannels=length(channels);

name0(1,1:7)=name(1,1:7);
name0(1,8:11)='conv';

name1(1,1:7)=name(1,1:7);
name1(1,8:11)='cnst';

fid0=fopen(name1,'r');
window=fread(fid0,1,'int32');
marg=fread(fid0,1,'int32');
fclose(fid0);

l_dane=raportconv(name0,header,nrchns); %dlugosc kanalu w *conv

if samples(1)<1
    error('Zbyt male wartosci w "samples"');
end

if samples(2)>l_dane
    error('Zbyt duze wartosci w "samples"');
end

%samples;
start=round((samples(1)-marg)/window);
stop=round((samples(2)-marg)/window)+1;

x=[((start-0.5)*window+marg):window:((stop-0.5)*window+marg)];

l_cnst=stop-start+1
y=zeros(lchannels,l_cnst);

lr=raportconv(name1,12,65);
start1=max(start,1);
start1=min(start1,lr-1);

stop1=min(stop,lr);
stop1=max(stop1,2);

l1=stop1-start1+1

y1=readconv(name1,12,nrchns,channels,[start1 stop1]);
size(y1)

if (start>=1)
    y(:,1:l1)=y1;
    y(:,l1:l_cnst)=y1(:,l1);
else
    y(:,(2-start):(2-start+l1-1))=y1;
    y(:,(2-start+l1):l_cnst)=y1(:,l1);
    y(:,1:(1-start))=y1(:,1);
end
clear y1;
size(x)
size(y)
c=interp1(x,y,[samples(1):samples(2)],'nearest');

s=readconv(name0,header,nrchns,channels,samples)-c;   