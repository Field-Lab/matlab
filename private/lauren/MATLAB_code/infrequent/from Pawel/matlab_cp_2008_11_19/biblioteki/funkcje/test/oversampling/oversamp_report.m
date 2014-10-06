function y=oversamp_report(signals,spikes,thresholds,filtr_dl,figura,osie,marg_left,marg_right);
%spikes: okresla NUMERY SPIKOW (sposrod znalezionych funkcja find!! - 
%patrz nizej) dla sygnalu s, ktore maja byc rysowane; thresholds -
%thresohld dla KAZDEGO spika
s_signals=size(signals);
dl=s_signals(1,2);
marg_left=15;
marg_right=20;
ts=[0:marg_right+marg_left-1]/20; %czas w milisekundach

figure(figura);
figura
for i=1:4
    i
    s=signals(i,:);
    
    y=blad_vs_energy(s,28); %szukamy spikow i energii bledu interpolacji
    a=find(y(3,:)>thresholds(1,i)); %sposrod znanych spikow wyszukujemy to o okreslonej energii
    la=length(a)
    thresholds(1,i)
    wskazniki=y(1,a) 
    
    signal=s(1:2:dl); %undersampling
    [y0,filtr]=oversampling(signal,2,filtr_dl,0.98); % oversampling;
    y1=y0(1,filtr_dl+1:filtr_dl+dl); %odrzucenie probek skrajnych
    roznica=y1-s;
    
    sp1_start=wskazniki(spikes(1,i))-marg_left;
    sp1_stop=wskazniki(spikes(1,i))+marg_right-1;
    spike_wspol=[sp1_start:sp1_stop];
    spike=s(spike_wspol);
    spike=spike-mean(spike);
    spike_y1=y1(spike_wspol)-mean(y1(spike_wspol));
    %qwer=mean(spike_y1)
    subplot(3,4,i);
    plot(ts,spike,'b-',ts,spike_y1,'r-');
    axis([0 1.7 -200 200]);
    if i==1        
        xlabel('time [ms]');
        ylabel('signal level');
    end
    if i==2                
        legend('original signal','interpolated signal');
    end
    grid on;
    
    subplot(3,4,i+4);
    plot(ts,spike,'bd-',ts,spike_y1,'rd-');
    axis([0.4 1 -200 200]);
    if i==2        
        xlabel('time [ms]');
        ylabel('signal level');
        legend('original signal','interpolated signal');
    end
    grid on;
    
    subplot(3,4,i+8);
    fs=fft(spike,500);
    f=[0:499]/500*20000;
    plot(f,abs(fs));
    xlabel('frequecy [Hz]');
    ylabel('Fourier coefficient');
    grid on;
    if i==1
        axis([0 10000 0 1000]);
    end
end