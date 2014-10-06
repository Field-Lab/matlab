czerwony_u=[0:0.01:0.33];
czerwony_i=[32:-1:28 26:-1:22 20:-1:17 15:-1:12 10 9 6.4 5.9 5.4 5 4.7 4.1 3.7 3.3 2.8 2.4 2 1.6 1.1 0.7 ];

n=21;
x=czerwony_i(1,1:n);
y=czerwony_u(1,1:n);

[p,s]=polyfit(x,y,1)
z=polyval(p,x);

plot(x,y,'bd',x,z);
%grid on;
e=1.6e-19;
h=6.63e-34;
l_cz=630e-9;
c=3e8;

E_cz=h*c/l_cz/e;
Ek=p(1,2);
W=E_cz-Ek;

zielony_i=[33:-1:21 19 18 16 14 12 11 9 8];
zielony_u=[0:0.01:0.12 0.14:0.02:0.28];

n=21;
x=zielony_i(1,1:n);
y=zielony_u(1,1:n);

[p,s]=polyfit(x,y,1)
z=polyval(p,x);

plot(x,y,'bd',x,z);
%grid on;
e=1.6e-19;
h=6.63e-34;
l_cz=630e-9;
c=3e8;

%E_cz=h*c/l_cz/e
Ek=p(1,2)
%W=E_cz-Ek
l_ziel=h*c/(Ek+W)/e




zielony_i=[77 70 62 54 47 41 35 29 24 19 14 12 10];
zielony_u=[0:0.04:0.36 0.4:0.02:0.44];

n=13;
x=zielony_i(1,1:n);
y=zielony_u(1,1:n);

[p,s]=polyfit(x,y,1)
z=polyval(p,x);

%plot(x,y,'bd',x,z);
axis([0 40 0 0.3]);
%grid on;
e=1.6e-19;
h=6.63e-34;
l_cz=630e-9;
c=3e8;

%E_cz=h*c/l_cz/e
Ek=p(1,2)
%W=E_cz-Ek
l_ziel=h*c/(Ek+W)/e



nieb_i=[51 46 43 39 35 31 27 24 20 17 16 14 13 12 11 9];
nieb_u=[0:0.04:0.36 0.38:0.02:0.48];

n=13;
x=nieb_i(1,1:n);
y=nieb_u(1,1:n);

[p,s]=polyfit(x,y,1)
z=polyval(p,x);

plot(x,y,'bd',x,z);
axis([0 60 0 0.5]);
%grid on;
e=1.6e-19;
h=6.63e-34;
l_cz=630e-9;
c=3e8;

%E_cz=h*c/l_cz/e
Ek=p(1,2)
%W=E_cz-Ek
l_ziel=h*c/(Ek+W)/e