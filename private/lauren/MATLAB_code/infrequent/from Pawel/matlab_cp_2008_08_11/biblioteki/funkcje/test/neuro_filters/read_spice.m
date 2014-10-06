function [ampl,phase]=read_spice(ampl_filename,phase_filename);
%Funkcja odczytuje dane z plikow wyeksportowanych ze SPICE, 
%zwraca dwie struktury dwukolumnowe (czest i ampl, czest i faza).
a=importdata(ampl_filename);
%spice podaje wartosci powyzej tysiaca troche dziwnie, np 1200
% to 1.2K. Matlab ignoruje 'K' i zotaje 1.2 Korekta ponizej:
%Dla czestosci = proste, musi monotonicznie rosnac.
czest1=a(:,1);
czest_d=diff(czest1);
gr=find(czest_d<0);
czest1(gr+1:length(czest1))=czest1(gr+1:length(czest1))*1000;
%Korekta wzmocnienia:
ampl=a(:,2);
ampl_d=diff(ampl);
gr1=find(ampl_d<-500)
size(gr1)
gr2=find(ampl_d>500)
if gr1
	if gr2
		ampl(gr1+1:gr2)=ampl(gr1+1:gr2)*1000;
	else
		ampl(gr1+1:length(ampl))=ampl(gr1+1:length(ampl))*1000;
	end
end
%plot(ampl)

a=importdata(phase_filename);
czest2=a(:,1);
czest_d=diff(czest2);
gr=find(czest_d<0);
czest2(gr+1:length(czest2))=czest2(gr+1:length(czest2))*1000;
phase=a(:,2);

ampl=[czest1 ampl];
phase=[czest2 phase];
