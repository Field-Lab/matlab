function [srd,odch]=channel_stat(dane);
%Funkcja zwraca dwa wektory, zawierajace odpowiednio wartosci srednia (srd)
%oraz odchylenia standardowe (odch) dla kolejnych kanalow tablicy dane. 
r=size(dane);

for i=1:r(1)
   srd(i)=mean(dane(i,:));
   odch(i)=std(dane(i,:)-srd(i));
end
