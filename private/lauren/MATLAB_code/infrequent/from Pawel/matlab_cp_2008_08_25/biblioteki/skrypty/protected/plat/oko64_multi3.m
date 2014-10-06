%oko64_multi.m
tic

k=0;

gain_aproks=1200; %orientacyjna wartosc, potrzebna dla podjecia decyzji czy w danym kanale 
%jest uzyteczny sygnal
%parametry poczatkowe:
ilosc=7;

cd /home/pawel/pliki/oko/oko_64/marzec2002
name(1)={'chipB_nast4_2czest2.dat'};
%name(2)={'chipB_gr2_all_czest.dat'};
%name(3)={'chipB_gr3_all_czest.dat'};
%name(4)={'chipB_gr4_all_czest.dat'};
%name(2)={'chipA_gr2'};
%name(3)={'chipA_gr3'};
%name(4)={'chipA_gr4'};
ec='wrgleirwhbgerg'
sname=size(name);

%dane=importdata('chipB_gr3_370Hz_2.dat');
%size(dane)

tlumik=2.12/0.026;
gen_ampl=[0.01:0.004:0.08];
ampl=0.02;
%ilosc=length(gen_ampl);
ilosc=20;

start=[0.1 0.05 -1 0.2];  %ampl,czest,faza,offset
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
%for plik=1:1
name(plik)
nazwa=name{plik}

dane=importdata(nazwa);

 for i=1:20
   nr_amplitudy=i
   %ampl=gen_ampl(i);
   
   prog=ampl/tlumik*gain_aproks*0.1;  %0.3 - wspolczynnik magiczny
   
   d=extract_multi(dane,i,64);
   
   ds=size(d);
   
   p=oko64_fit(d,'sinus_fit',start,krok,prog);  %   
   
   'ccc'
   
   w1=p(:,1);
   w2=p(:,2);
   %size(w);
   e=find((abs(w1)>prog) & (abs(w2)>0.025) & (abs(w2)<0.1))  %detekcja kanalow z sygnalem
  % channels=[channels e'];
   
   gain(e,i)=p(e,1);    %wpisanie wzmocnienia dla danej amplitudy wejsciowej i dla tych...
   e;
                        %...kanalow, ktore uznano za aktywne
   offset(e,i)=p(e,4);  %to co powyzej, dla offsetow
   e
   
   for j=1:length(e)
      j
      kanal=e(j,1)
      d(kanal,:);
      y=feval('sinus_fit',p(kanal,1:4),[1:ds(2)]);
      ys=size(y);
      q=wsp_kor(d(kanal,:),y);
      korr(kanal,i)=q;
      %min(d(kanal,:))
      %min(y)
      %max(d(kanal,:))
      %max(y)
      if q<0.9999
        i
        k=k+1;
        wpol_korelacji=q;
        %min_signal=min(d(e(j,1),:))
	%min_fit=min(y)
	%max_signal=max(d(e(j,1),:))
	%max_fit=max(y)
	%size(d(kanal,:))
	figure(i)
	subplot(8,8,j)
	plot([1:32],d(kanal,:),[1:32],y);
	warning('Mala wartosc wspolczynnika korelacji');
      end
   end
 end % of i
 
end % of plik
%channels=[channels e'];

n1=name{1};
n2=n1(1:5);

n3=[n2 '_all_gain_nast4.dat'];
save(n3,'gain','-ASCII','-DOUBLE');

n3=[n2 '_all_offset_nast4.dat'];
save(n3,'offset','-ASCII','-DOUBLE');

n3=[n2 '_all_korr_nast4.dat'];
save(n3,'korr','-ASCII','-DOUBLE');

setup=[tlumik gen_ampl];
n3=[n2 '_all_setup_nast4.dat'];
save(n3,'setup','-ASCII','-DOUBLE');

toc
