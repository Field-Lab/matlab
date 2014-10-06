function [h1] = filtracja (t,x,fcut)
% t - wektor czasu
% s - sygna³
% fcut - czêstotliwoœæ graniczna

%  wzorowane na kodzie: Shmuel Ben-Ezra, Ultrashape ltd. August 2009

%% weryfikacja danych
if ~any(size(t)==1),
   disp('Z³y format!Powinna byæ macierz jednowymiarowa.')
   return
end
if ~any(size(x)==1),
   disp('Z³y format!Powinna byæ macierz jednowymiarowa.')
   return
end
if length(t)~=length(x),
   disp('Z³y foramt! Wektory powinny mieæ te same d³ugoœci.')
   return
end
%% Definicje
Fs=1/(t(2)-t(1)); %czêstotliwoœæ próbkowania
N=length(x);
Nfft=2^nextpow2(N);
f=Fs/2*linspace(0,1,1+Nfft/2); % wektor czêstotliwoœæi
cutoff_freq=fcut; 
my_freqs=[];
%% g³ówny kod
y=fft(x,Nfft)/N; % transformata 
y2=filterfft(f, y, cutoff_freq, my_freqs); % filtrowane czêst.
X=ifft(y2); % odwrotna transformata
X=X(1:N)/max(X);
ind1 = find(y2(1:1+Nfft/2)); % znalezienie niezerowych elementów
nf1 = length(ind1); % liczenie niezerowych elementów
%% wykresy
figname = 'FFT';
ifig = findobj('type', 'figure', 'name', figname);
if isempty(ifig),
    ifig = figure('name', figname); 
end
figure(ifig);
% pierwszy wykres
subplot(3,1,1)
plot(t,x)
xlabel('Sec')
axis tight
title('Oryginalny sygna³')
%drugi wykres
subplot(3,1,2)
yplot=abs(y(1:1+Nfft/2));
yplot=yplot/max(yplot);
semilogy(f*1e-3, yplot, f(ind1)*1e-3, yplot(ind1), '.r');
xlabel('kHz')
title('Amplitudy')
legend('ca³e widmo', 'wybrane czêstotliwoœci')
% trzeci wykres
subplot(3,1,3)
plot(t,X)
xlabel('Sec')
if isempty(cutoff_freq),
    scutoff='Brak czêstotliwoœci granicznej';
else
    scutoff=sprintf('Czêst. graniczna = %g [khz]', cutoff_freq/1e3);
end
stitle3=sprintf('Zrekonstruowany sygna³ z %d wybranych czêstotliwoœci; %s', nf1, scutoff);
title(stitle3)
axis tight
return

function y2=filterfft(f, y, fcut, wins)
nf=length(f);
ny=length(y);
if ~(ny/2+1 == nf),
    disp('Z³e wymiary wprowadzonych danych!')
    y2=-1;
    return
end

% filtr czêstotliwoœci
y2=zeros(1,ny);
if ~isempty(fcut)
    ind1=find(f>=fcut);
    y2(ind1) = y(ind1); % 
else
    y2=y;
end

return