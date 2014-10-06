start=1000000;
dlugosc=2500000;

%name='Data013';
channel=65;

prog=400;
histereza=50;

figura=13;

cd /mnt/win3/oko/2000-12-12;

figure(figura);
clf(figura);

ax=[5 100000 300 1200];
st1=1000000;
st2=6000000;

read_param=struct('name','Data009','header',206,'nrchns',65,'samples',[start+1 start+dlugosc]);


znak=1;
%plot(odl,wys2,'bd');
name='Data009';
start=st1;
[odl,wys]=ampl_vs_opozn(name,206,65,channel,[start+1 start+dlugosc],prog,histereza,znak);
subplot(2,3,1);
if length(odl)>1
    semilogx(odl,wys,'bd');
end
axis([ax(1) ax(2) ax(3) ax(4)]);
grid on;

start=st2;
[odl,wys]=ampl_vs_opozn(name,206,65,channel,[start+1 start+dlugosc],prog,histereza,znak);
subplot(2,3,4);
if length(odl)>1
    semilogx(odl,wys,'bd');
end
axis([ax(1) ax(2) ax(3) ax(4)]);
grid on;

name='Data011';
start=st1;
[odl,wys]=ampl_vs_opozn(name,206,65,channel,[start+1 start+dlugosc],prog,histereza,znak);
subplot(2,3,2);
if length(odl)>1
    semilogx(odl,wys,'bd');
end
axis([ax(1) ax(2) ax(3) ax(4)]);
grid on;

start=st2;
[odl,wys]=ampl_vs_opozn(name,206,65,channel,[start+1 start+dlugosc],prog,histereza,znak);
subplot(2,3,5);
if length(odl)>1
    semilogx(odl,wys,'bd');
end
axis([ax(1) ax(2) ax(3) ax(4)]);
grid on;

name='Data013';
start=st1;
[odl,wys]=ampl_vs_opozn(name,206,65,channel,[start+1 start+dlugosc],prog,histereza,znak);
subplot(2,3,3);
if length(odl)>1
    semilogx(odl,wys,'bd');
end
axis([ax(1) ax(2) ax(3) ax(4)]);
grid on;

start=st2;
[odl,wys]=ampl_vs_opozn(name,206,65,channel,[start+1 start+dlugosc],prog,histereza,znak);
subplot(2,3,6);
if length(odl)>1
    semilogx(odl,wys,'bd');
end
axis([ax(1) ax(2) ax(3) ax(4)]);
grid on;
'koniec.'