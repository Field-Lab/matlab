f1=8996.4;
A=0.05;

fs=44100;
N=16;
ts=1/(N*fs);

n=1000;

xi=0;
ai=0;
bi=0;
yi=0;

len=N*n;
y=zeros(1,len);

for i=1:len
    i;
    xi=A*sin(2*pi*f1*i*ts);
    ai=xi-yi;
    bi=bi+ai;
    y0=sign(bi);
    if y0~=0
        yi=y0;
    end
    y(i)=yi;
end

subplot(2,1,1);
plot(y);
fy=fft(y)*2/len;
f=[0:len-1]/len*N*fs;
subplot(2,1,2);
plot(f,abs(fy));
axis([0,N*fs/2,0,1]);
grid on;