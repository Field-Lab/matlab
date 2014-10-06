%oko64_multi.m
tic

k=0;

gain_aproks=500; %orientacyjna wartosc, potrzebna dla podjecia decyzji czy w danym kanale 
%jest uzyteczny sygnal
%parametry poczatkowe:
ilosc=7;

cd /home/pawel/pliki/oko/oko_64/
dane=importdata('chipB_gr3_370Hz_2.dat');
size(dane)

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

for i=1:length(gen_ampl)
   nr_amplitudy=i
   ampl=gen_ampl(i);
   prog=ampl/tlumik*gain_aproks*0.3;
   
   d=extract_multi(dane,i,64);
   ds=size(d);
   p=oko64_fit(d,'sinus_fit',start,krok,prog);  %
   
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
	figure(i)
	subplot(3,3,j)
	plot([1:32],d(kanal,:),[1:32],y);
	warning('Mala wartosc wspolczynnika korelacji');
      end
   end
      
     
   %gain=[gain w(e)'];
   %w=p(:,4);
   %offset=[offset w(e)'];
     
end
%channels=[channels e'];



%kolory=['b','k','y','m','c','r','g'];

%figure(1);
%for i=1:ilosc
%  plot(gen_ampl,gain(:,i),[kolory(i) 'd-']);
%  hold on;
%end

%figure(2);
%for i=1:ilosc
%  plot(gen_ampl,offset(:,i),[kolory(i) 'd-']);
%  hold on;
%end


toc
