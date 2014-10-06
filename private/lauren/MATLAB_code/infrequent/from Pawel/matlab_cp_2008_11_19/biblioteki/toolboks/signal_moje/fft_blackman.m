function [s_tb,f]=fft_blackman(s,marg,N);
%[s_tb,f]=fft_blackman(s,marg,N);
%Funkcja zwraca widmo fft sygnalu przemnozonego przez funkcje okna
%trapez_blackman. Argumenty: s - sygnal, marg - margines (po lewej i
%prawej) definiujacy przedzialy probek mnozonych przez wartosci mniejsze 0d
%1; N - ilosc zwracanych wspolczynnikow fft (jesli N==0, rowna
%ilosci probek);
%Funkcja zwraca wartosc sygnalu po "zokienkowaniu" (s_tb) oraz widmo.
tb=trapez_blackman(length(s),marg);
s_tb=s.*tb;
if N==0
    f=fft(s_tb);
else
    f=fft(s_tb,N);
end