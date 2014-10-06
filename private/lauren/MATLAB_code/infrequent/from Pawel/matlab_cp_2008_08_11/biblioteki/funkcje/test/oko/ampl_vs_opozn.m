function [odl,wys]=ampl_vs_opozn(read_param,detect_param);
%function [odl,wys]=ampl_vs_opozn(name,header,nrchns,channel,samples,prog,histereza);
%Funkcja dokonujaca zestawienia: apmlituda piku vs jego opoznienie
%w czasie w stosunku do piku poprzedniego. 

w=detect2_f(read_param,detect_param);

ilosc=length(w);

%'analyzing...'
if ilosc>1
   pocz=w(1,:);
   kon=w(2,:);
   wys=w(3,:);

   odl=diff(pocz);
   szer=kon(1,2:length(w))-pocz(1,2:length(w));
   wys=abs(wys(1,2:length(w)));
else
    odl=0;
    wys=0;
end
