fs = 2000;             %cz�stotliwo�� pr�bkowania
f0 = 60;                %cz�stotliwo�� notch (do usuni�cia)
fn = fs/2;              %cz�stotliwo�� Nyquista
freqRatio = f0/fn;      
t=0:0.0005:1;
notchWidth = 0.01;       %szeroko�� notch
L=length(t);
%Policzenie zer
zeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

%Policzenie biegun�w
poles = (1-notchWidth) * zeros;

%figure;
%zplane(zeros.', poles.');

b = poly( zeros ); 
a = poly( poles );

%figure;
%freqz(b,a,32000,fs)
x=3*sin(2*pi*60*t)+3*sin(2*pi*100*t)+3*sin(2*pi*50*t);
%filtracja sygna�u x
y = filter(b,a,x);



 m = length(x);          % d�ugo�� sygna�u
n = pow2(nextpow2(m));  % d�ugo�� transformaty
Y = fft(x,n);           % DFT
f = (0:n-1)*(fs/n);     % wektor cz�stotliwo��i
power = Y.*conj(Y)/n;   % Moc DFT


Y2 = fft(y,n);           % DFT
power2 = Y2.*conj(Y2)/n;   % Moc DFT

figure;
subplot(4,1,1);
plot(t,x)
subplot(4,1,2)
plot(f,power) 
title('Moc DFT x')
xlabel('Cz�stotliwo�� (Hz)')
ylabel('Moc')
subplot (4,1,3);
plot(t,y)
subplot (4,1,4);
plot(f,power2) 
title('Moc DFT y')
xlabel('Cz�stotliwo�� (Hz)')
ylabel('Moc')
