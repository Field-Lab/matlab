%oko64_multi.m
tic

k=0;

gain_aproks=500; %orientacyjna wartosc, potrzebna dla podjecia decyzji czy w danym kanale 
%jest uzyteczny sygnal
%parametry poczatkowe:
ilosc=7;

cd /home/pawel/pliki/oko/oko_64/koncowe
name(1)={'chipA_gr1_czest_2.dat'};
name(2)={'chipA_gr2_czest_2.dat'};
name(3)={'chipA_gr3_czest_2.dat'};
name(4)={'chipA_gr4_czest_2.dat'};
%name(2)={'chipA_gr2'};
%name(3)={'chipA_gr3'};
%name(4)={'chipA_gr4'};
sname=size(name);

%dane=importdata('chipB_gr3_370Hz_2.dat');
%size(dane)

tlumik=2.12/0.026;
gen_ampl=[0.01:0.004:0.08];
ilosc=length(gen_ampl);

start=[0.1 0.05 -1 0.2];
krok=[0.0001 0.0002 0.0005 0.0002];
prog=0.05;

%gain=zeros(ilosc,10);
%offset=zeros(ilosc,10);
gain=zeros(64,ilosc);
f=gain;
faza=gain;
offset=gain;
korr=gain;
channels=[];

for plik=1:sname(2)
name(plik)
'ksdfgkjehfgvg'
nazwa=name{plik}
dane=importdata(nazwa);
ec='jjjjj'
 for i=1:length(gen_ampl)
   nr_amplitudy=i
   ampl=gen_ampl(i);
   prog=ampl/tlumik*gain_aproks*0.3;  %0.3 - wspolczynnik magiczny
   
   'rhrtyj56yj6u5jk'
   
   i
   d=extract_multi(dane,i,64);
   prog=max((max(d')-min(d'))/2*0.4)
   
   'kjkjkjkjkjk'
   
   ds=size(d);
   p=oko64_fit(d,'sinus_fit',start,krok,prog);  %
   
   'hhhhhhhhh'
   
   w=p(:,1);
   size(w);
   e=find(abs(w)>prog);  %detekcja kanalow z sygnalem
  % channels=[channels e'];
   
   gain(e,i)=p(e,1);    %wpisanie wzmocnienia dla danej amplitudy wejsciowej i dla tych...
   e;
                        %...kanalow, ktore uznano za aktywne
   offset(e,i)=p(e,4);  %to co powyzej, dla offsetow
   e;
   
   for j=1:length(e)
      kanal=e(j,1)
      d(kanal,:)
      y=feval('sinus_fit',p(kanal,1:4),[1:ds(2)]);
      ys=size(y);
      q=wsp_kor(d(kanal,:),y);
      korr(kanal,i)=q;
      %min(d(kanal,:))
      %min(y)
      %max(d(kanal,:))
      %max(y)
      if q<0.9995
        i
        k=k+1;
        wpol_korelaji=q
        %min_signal=min(d(e(j,1),:))
	%min_fit=min(y)
	%max_signal=max(d(e(j,1),:))
	%max_fit=max(y)
	%size(d(kanal,:))
	figure(i+plik*20)
	subplot(3,3,j)
	plot([1:32],d(kanal,:),[1:32],y);
	warning('Mala wartosc wspolczynnika korelacji');
      end
   end
 end % of i
 
end % of plik
%channels=[channels e'];

n1=name{1};
n2=n1(1:5);

n3=[n2 '_gain.dat'];
save(n3,'gain','-ASCII','-DOUBLE');

n3=[n2 '_offset.dat'];
save(n3,'offset','-ASCII','-DOUBLE');

n3=[n2 '_korr.dat'];
save(n3,'korr','-ASCII','-DOUBLE');

setup=[tlumik gen_ampl];
n3=[n2 '_setup.dat'];
save(n3,'setup','-ASCII','-DOUBLE');

toc
