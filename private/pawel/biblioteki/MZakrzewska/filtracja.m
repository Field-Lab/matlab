function [h1] = filtracja (t,x,fcut)
% t - wektor czasu
% s - sygna�
% fcut - cz�stotliwo�� graniczna

%  wzorowane na kodzie: Shmuel Ben-Ezra, Ultrashape ltd. August 2009

%% weryfikacja danych
if ~any(size(t)==1),
   disp('Z�y format!Powinna by� macierz jednowymiarowa.')
   return
end
if ~any(size(x)==1),
   disp('Z�y format!Powinna by� macierz jednowymiarowa.')
   return
end
if length(t)~=length(x),
   disp('Z�y foramt! Wektory powinny mie� te same d�ugo�ci.')
   return
end
%% Definicje
Fs=1/(t(2)-t(1)); %cz�stotliwo�� pr�bkowania
N=length(x);
Nfft=2^nextpow2(N);
f=Fs/2*linspace(0,1,1+Nfft/2); % wektor cz�stotliwo��i
cutoff_freq=fcut; 
my_freqs=[];
%% g��wny kod
y=fft(x,Nfft)/N; % transformata 
y2=filterfft(f, y, cutoff_freq, my_freqs); % filtrowane cz�st.
X=ifft(y2); % odwrotna transformata
X=X(1:N)/max(X);
ind1 = find(y2(1:1+Nfft/2)); % znalezienie niezerowych element�w
nf1 = length(ind1); % liczenie niezerowych element�w
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
title('Oryginalny sygna�')
%drugi wykres
subplot(3,1,2)
yplot=abs(y(1:1+Nfft/2));
yplot=yplot/max(yplot);
semilogy(f*1e-3, yplot, f(ind1)*1e-3, yplot(ind1), '.r');
xlabel('kHz')
title('Amplitudy')
legend('ca�e widmo', 'wybrane cz�stotliwo�ci')
% trzeci wykres
subplot(3,1,3)
plot(t,X)
xlabel('Sec')
if isempty(cutoff_freq),
    scutoff='Brak cz�stotliwo�ci granicznej';
else
    scutoff=sprintf('Cz�st. graniczna = %g [khz]', cutoff_freq/1e3);
end
stitle3=sprintf('Zrekonstruowany sygna� z %d wybranych cz�stotliwo�ci; %s', nf1, scutoff);
title(stitle3)
axis tight
return

function y2=filterfft(f, y, fcut, wins)
nf=length(f);
ny=length(y);
if ~(ny/2+1 == nf),
    disp('Z�e wymiary wprowadzonych danych!')
    y2=-1;
    return
end

% filtr cz�stotliwo�ci
y2=zeros(1,ny);
if ~isempty(fcut)
    ind1=find(f>=fcut);
    y2(ind1) = y(ind1); % 
else
    y2=y;
end

return