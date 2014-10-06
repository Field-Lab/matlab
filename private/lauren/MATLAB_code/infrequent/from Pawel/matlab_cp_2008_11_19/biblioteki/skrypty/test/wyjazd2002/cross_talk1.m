frqs=[330 700 1000 1600 2000 2700];
channel_512=[1260 1280 1240 1120 1040 880];
channel_511=[4 7 7 8 9 9.5];

a=channel_511./channel_512

a=20*log10(a)

figure(2);
semilogx(frqs,a,'k-*');
grid on;