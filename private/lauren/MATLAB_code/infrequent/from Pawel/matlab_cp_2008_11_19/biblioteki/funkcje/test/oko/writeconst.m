function y=writeconst(name,window,marg,prog,header,nrchns);

name0=name(1,1:7)
name1(1,1:7)=name0;
name1(1,8:11)='conv'
fid0 = fopen(name1,'r');

namesize=size(name);
fileoutname(1,1:7)=name0;
fileoutname(1,8:11)='cnst'
fid1=fopen(fileoutname,'w');

a=fseek(fid0,0,1); %skok na koniec pliku
p=ftell(fid0);
size=p;

fclose(fid0);

chlength=floor((p-header)/2/nrchns); %samples per channel

fwrite(fid1,window,'int32');
fwrite(fid1,marg,'int32');
fwrite(fid1,prog,'int32');
header=12;

for i=1:nrchns
    i
    r3=readconv(name1,header,nrchns,[i],[1 chlength]);
    %length(r3)
    %size(r3)
    b=round(shwconst2(r3,window,marg,prog));
    %p=ftell(fid1)
    fwrite(fid1,b','int16');
    l=ftell(fid1);
end

fclose(fid1);
y=length(b);