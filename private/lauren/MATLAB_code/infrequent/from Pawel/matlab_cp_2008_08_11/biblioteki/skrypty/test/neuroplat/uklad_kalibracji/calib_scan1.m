%Pomiar kalibracji - skanowanie DACa oraz kanalow dla DAC=15. Kanal nr 2 (numeracja od zera).

%1.Liczenie
calDAC=[15.5 30.5 46 61 76 91.6 106 120 136 150 165 181 194 209 224 238 253 268 282 296 311 326 339 354 367 382 396 410 426 440 454]*10; %in uV
calChns=[223 223 224 224 225 225 225 225 226 226 226 226 227 227 227 227 227 227 227 227 228 228 228 228 228 228 227 228 228 228 229 229 229 228 228 228 228 228 228 227 227 227 227 227 227 226 226 226 226 226 225 225 224 224 224 224 223 223 222 222 221 221 220 220]*10; %in uV

mnoznik=3.28/1.43*1e-3;
calDAC=calDAC*mnoznik;
calChns=calChns*mnoznik;

%2.Rysowanie
textsize=18;

figure(1);
plot(calDAC,'bd-');
grid on;
axis([0 32 0 12]);
h1=gca;
set(h1,'FontSize',textsize);
h1=xlabel('calib DAC value [LSB]');
set(h1,'FontSize',textsize);
h1=ylabel('voltage P-P [mV]')
set(h1,'FontSize',textsize);

figure(2);
plot(calChns,'bd-');
grid on;
axis([0 65 5 5.3]);
h1=gca;
set(h1,'FontSize',textsize);


h1=xlabel('channel number');
set(h1,'FontSize',textsize);
h1=ylabel('voltage P-P for DAC=15 [mV]');
set(h1,'FontSize',textsize);
tx=52;
ty=5.28;
h1=text(tx,ty,['mean:' num2str(mean(calChns),3)]);
set(h1,'FontSize',textsize);
h1=text(tx,ty-0.02,['std/mean:' num2str(std(calChns)/mean(calChns)*100,'%4.2f') '%']);
set(h1,'FontSize',textsize);
