cd /home/pawel/pliki/oko/oko_64/15-01-2002/

tlumik=2.12/0.026;
gen_ampl=0.02;

name(1)={'chipA_gr1'};
name(2)={'chipA_gr2'};
name(3)={'chipA_gr3'};
name(4)={'chipA_gr4'};
s_name=size(name)
start=[0.1 0.05 0.0 0.2];
krok=[0.0001 0.0002 0.0005 0.0002];

ch_tekst='channel: ';
gain_tekst='gain: ';
offset_tekst='offset: ';

prog=0.05;;
kanaly=zeros(0);
przebiegi=zeros(0);
parametry=zeros(0);

for i=1:s_name(2)
  name{i}
  q=oko64_extract(name{i},64);
  p=oko64_fit(q,'sinus_fit',start,krok,prog);  %
  %p(:,1)=p(:,1)*tlumik/gen_ampl);
  %wyznaczenie tych kanalow, gdzie byl podawany sygnal
  w=p(:,1);
  e=find(abs(w)>prog);
  for j=1:length(e)
    przebiegi=[przebiegi' q(e(j),:)']';
    parametry=[parametry' p(e(j),:)']';
    size(przebiegi)
    size(parametry)
  end
  kanaly=[kanaly e'];
end

figure(2);
sq=size(q);
t=[1:sq(2)];
clear sq;
%kanaly=sort(kanaly);
n0=ceil(sqrt(length(kanaly)));
n1=ceil(length(kanaly)/n0);

for i=1:length(kanaly)
  subplot (n0,n1,i);
  s=przebiegi(i,:);
  a=feval('sinus_fit',parametry(i,:),t);
  plot(t,s,t,a);
  axis([1 32 0 0.5]);
  
  channel=num2str(kanaly(i));
  text(10,0.45,[ch_tekst channel]);
  
  ampl=num2str(abs(parametry(i,1))*tlumik/gen_ampl,'%2.3f')
  text(10,0.38,[gain_tekst ampl]);
  
  offset=num2str(parametry(i,4),'%2.3f')
  text(10,0.31,[offset_tekst offset]);
  grid on;
end

%stare=[577 585 445 630 645 649 570 640 638 570 634 562 619 622 576 

