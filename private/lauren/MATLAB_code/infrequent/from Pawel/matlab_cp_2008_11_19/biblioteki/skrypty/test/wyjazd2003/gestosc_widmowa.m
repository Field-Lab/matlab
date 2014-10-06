cd /home2/pawel/oko/santa_cruz_2003/2003_10_21;
filename='nostim_noise_water_Vpol-14Vconv';
fs=20000;
start=1;
dlugosc=100*fs;
lsb=1;

channels=[1:4];

figure(4);
for i=1:4
    s=readconv(filename,206,65,channels(i),[(start) (start+dlugosc-1)]);
    s=s-mean(s);
    [f,w]=spdf(s,20000,20000,fs,lsb);
    subplot(2,2,i);
    loglog(w);
    axis([1 round(fs/2) 1e-4 1e4]);
    grid on;
    suma=sum(w)
    sig=sqrt(suma/2);
    text(1000,1000,num2str(sig));
    suma2=suma-w(61)-w(181)-w(301);
    sig2=sqrt(suma2/2);
    text(1000,300,num2str(sig2));
end
