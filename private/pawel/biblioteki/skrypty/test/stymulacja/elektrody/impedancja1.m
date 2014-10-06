n=10; %ilosc elektrod
m=19; %ilosc czest.

f=zeros(n,m);
z=f;
fi=f;

% * * * pomiary 2005_09_20

f1=[20 30 40 70 130 270 380 520 670 1020 1320 1720 2020 2520 3020 4020 5020 7000 10000];
z1=[1170 944 818 610 440 290 223 183 156 120 103 88.5 81 71 64 55 50 49 45];
fi1=[-46 -46 -45 -50 -48 -47 -51 -51.5 -52 -51.5 -50.5 -50 -48.5 -48.5 -48 -47 -44 36 -32];


f2=[20 30 40 70 130 270 380 520 670 1020 1320 1720 2020 2520 3020 4020 5020 7000 10000];
z2=[2060 1700 1500 1150 830 557 454 373 319 245 202 170 156 135 119 98 85 68.5 55];
fi2=[-40.5 -42 -40 -46 -48 -52 -53 -54.5 -55.5 -56 -58 -58 -56.5 -56.5 -56 -56 -55.5 -54 -51.5];

% * * * pomiary 2005_09_21

f3=[20 30 40 60 90 130 270 520 1020 2020 4000 7000 10000 20000 40000 80000 150000];
z3=[1500 1210 1050 835 582 479 290 192 120 89 59 43.5 36 25 19 15 14];
fi3=-[48.1 49.5 52.2 52.9 56.6 56.8 57.6 55 54.2 55 51.5 49 46.5 42.4 36 31.5 27];

f4=f3;
z4=[3300 2440 2070 1600 1210 970 618 393 224 160 104 73 58 37 23.5 16 11.7];
fi4=-[57.2 55.5 56.4 55.8 62.5 58.5 60 60 60.5 58.6 58 57.6 57.5 56 53 46.9 41.4];

loglog(f1,z1,'bd-',f2,z2,'gd-',f3,z3,'kd-',f4,z4,'rd-');
grid on;
