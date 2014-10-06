header=108;
frame=770; %ramka - dane dla 512 kanalow + 2 bajty ekstra

cd I:\crosstalk_Sasha\data024;

fid = fopen('data024000.bin','r');
fseek(fid,header+2,-1);

N=2^10;
mnoznik=1.5;
F=zeros(512*mnoznik,N);

for i=1:N
    %a=ftell(fid)
    F(:,i)=fread(fid,512*mnoznik,'ubit8',0);
    fseek(fid,2,0);
    %a=ftell(fid)
end
fclose(fid);

offset=382*mnoznik;
b1=F(offset+1,:);
b2=F(offset+2,:);
b3=F(offset+3,:);

s1=b1*16+floor(b2/16)-2048;
s2=(b2-floor(b2/16)*16)*256+b3-2048;

fs=zeros(512,N);
plot(s2)

for i=1:256
    i
    offset=2*mnoznik*(i-1);
    b1=F(offset+1,:);
    b2=F(offset+2,:);
    b3=F(offset+3,:);

    s1=b1*16+floor(b2/16)-2048;
    s2=(b2-floor(b2/16)*16)*256+b3-2048;
    
    fs(2*i-1,:)=fft(s1);
    fs(2*i,:)=fft(s2);
end

figure(3)
semilogy(abs(fs(:,11))/abs(fs(384,11))*100);
grid on;