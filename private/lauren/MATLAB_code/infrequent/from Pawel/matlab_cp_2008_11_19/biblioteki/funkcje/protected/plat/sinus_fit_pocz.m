function parametry=sinus_fit_pocz(s);
%lot(s)

ampl=(max(s)-min(s))/2;
ofst=(max(s)+min(s))/2;
s=s-ofst;
figure(55)
plot(s)
ss=sign(s);
dss=diff(ss);

a=find(dss~=0); % gdzie zmienia sie znak funkcji - przejscia przez zero (offset)
if length(a)<2
   parametry=[0 0 0 ofst];
   dt=length(s);
else
  da=diff(a);  % odleglosci miedzy kolejnymi miejscami zerowymi
  mda=max(da); % szacunkowa oldleglosc na osi czasu pomiedzy dwoma kolejnymi
             % miejscami zerowymi sinusa - a nie tymi wynikajacymi z szumow
  prog=mda/2;   % 
  pts=find(da>prog);  % te odleglosci miedzy punktami zerowymi, ktore
                      % odpowiadaja faktycznie nastepnemu miejscu zerowemu 
  pts=[0 pts length(a)];
  DT=zeros(1,length(pts)-1);
  for i=1:length(pts)-1
    DT(i)=mean(a((pts(i)+1):pts(i+1)));
  end

  %DT
  dt=mean(diff(DT));
  %dt=a(1,2)-a(1,1)
  %'tutaj...'
  dt;
  faproks=1/(2*dt);
  %'czy tu'
  size(a);
  size(s);

  znak=sign(s(max(a(1,1)-1,1)));
  %'czy tez tu'
  size(DT);
  faza=(dt-DT(1))/dt*pi*znak;
  %'a moze tu'
  ofst;
  parametry=[ampl faproks faza ofst];
end
