%oko64_multi.m
tic

gain_aproks=500; %orientacyjna wartosc, potrzebna dla podjecia decyzji czy w danym kanale 
%jest uzyteczny sygnal
%parametry poczatkowe:
ilosc=7;

cd /home/pawel/pliki/oko/oko_64/
dane=importdata('chipB_gr1_370Hz_2.dat');
size(dane)

tlumik=2.12/0.026;
gen_ampl=[0.01:0.004:0.03];


start=[0.1 0.05 -1 0.2];
krok=[0.0001 0.0002 0.0005 0.0002];
prog=0.05;

gain=zeros(ilosc,10);
offset=zeros(ilosc,10);
gain=[];
f=[];
faza=[];
offset=[];
korr=[];


for i=1:length(gen_ampl)
   amplituda=i
   prog=3*gen_ampl(1,i);
   d=extract_multi(dane,i,64);
   %size(d)
   p=oko64_fit(d,'sinus_fit',start,krok,prog);  %
   w=p(:,1);
   e=find(abs(w)>prog);
   gain(i,1:length(e))=w(e)';
   w=p(:,4);
   %e=find(abs(w)>prog);
   offset(i,1:length(e))=w(e)';
end

kolory=['b','k','y','m','c','r','g'];

figure(1);
for i=1:ilosc
  plot(gen_ampl,gain(:,i),[kolory(i) 'd-']);
  hold on;
end

figure(2);
for i=1:ilosc
  plot(gen_ampl,offset(:,i),[kolory(i) 'd-']);
  hold on;
end


toc
