header=108;
frame=770; %ramka - dane dla 512 kanalow + 2 bajty ekstra

cd I:\crosstalk_Sasha\data037;

fid = fopen('data037000.bin','r');
fseek(fid,header+2,-1);

N=1024;
F=zeros(512,N);

for i=1:N
    %a=ftell(fid)
    F(:,i)=fread(fid,512,'bit12',0);
    fseek(fid,2,0);
    %a=ftell(fid)
end
fclose(fid);

channel=385;
s=F(channel,:);
s0=floor(s/256);
s1=s-s0*256;

figure(1);
subplot(4,2,1);
plot(s0,'bd-');
subplot(4,2,2);
plot(abs(fft(s0-mean(s0))));
subplot(4,2,3);
plot(s1,'bd-');
subplot(4,2,4);
plot(abs(fft(s1-mean(s1))));
%subplot(4,1,3);

s=F(channel+1,:);
s0=floor(s/256);
s1=s-s0*256;
subplot(4,2,5);
plot(s0,'bd-');
subplot(4,2,6);
plot(abs(fft(s0-mean(s0))));
subplot(4,2,7);
plot(s1);
subplot(4,2,8);
plot(abs(fft(s1-mean(s1))));


