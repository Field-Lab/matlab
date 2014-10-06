function y=crosstalk1(filename,fs,f0,channel0,wykladnik);

length=600000;
mnoznik=1.5;
header=108;

N=2^wykladnik;
steps=floor(length/N)
f_nr=round(f0/fs)*N+1; % nr czestotliwosci f0 w tr. fouriera

fid = fopen(filename,'r');
fseek(fid,header+2,-1);

F=zeros(512*mnoznik,N);

y=zeros(1,512);
for i=1:steps
    i
    %czytanie danych - N probek
    for j=1:N
        F(:,j)=fread(fid,512*mnoznik,'ubit8',0);
        fseek(fid,2,0);
    end
    
    %dekodowanie kanalow (N probek) i fft
    for i=1:256    
        offset=2*mnoznik*(i-1);
        b1=F(offset+1,:);
        b2=F(offset+2,:);
        b3=F(offset+3,:);

        s1=b1*16+floor(b2/16)-2048;
        s2=(b2-floor(b2/16)*16)*256+b3-2048;
    
        fs(2*i-1,:)=fft(s1);
        fs(2*i,:)=fft(s2);
    end
    %korekta fazy
    wspol=fs(channel0,f_nr);
    korekta=conj(wspol)/abs(wspol);
    dane=fs(:,nr0)*korekta;
    y=y+abs(dane');
end