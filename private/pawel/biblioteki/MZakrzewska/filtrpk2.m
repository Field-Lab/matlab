F=0:1/10000:1;
H2=[zeros(1,201),ones(1,9800)];
h2=fir2(200,F,H2);
figure
freqz(h2,1)