H=[zeros(1,201),ones(1,19599),zeros(1,200)];
h=real(ifft(H));
w=hamming(201)';
w20000=[w(101:201),zeros(1,19799),w(1:100)];
hw=h.*w20000;
h1=[hw(1:100),hw(19900:20000)];
h1=fftshift(h1);
freqz(h1,1);