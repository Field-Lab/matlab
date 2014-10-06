function y=oversamp_report_onechannel(s,detect_results,spikes,filter_param,figura,osie,margins,RC_filter);
%Robi analize dla kilku spikow, spikes: wynik funkcji detect (poczatki
%okreslonych spikow).

dl=length(s);
marg_left=margins(1,1);;
marg_right=margins(1,2);
ts=[0:marg_right+marg_left-1]/20; %czas w milisekundach
ls=length(spikes);

signal=s(1:2:dl); %undersampling
[y0,filtr]=oversampling(signal,filter_param.N,filter_param.order,filter_param.freq); % oversampling;
y1=y0(1,filter_param.order+1:filter_param.order+dl); %odrzucenie probek skrajnych
clear y0;
roznica=y1-s;

y2=conv(s,RC_filter);
lr=(length(RC_filter)-1)/2;
y2=y2(1,1+lr:dl+lr);
size(y2)

figure(figura);
for i=1:ls         
    sp1_start=detect_results(1,spikes(1,i))-marg_left
    sp1_stop=detect_results(spikes(1,i))+marg_right-1
    spike_wspol=[sp1_start:sp1_stop];
    spike=s(spike_wspol);
    spike=spike-mean(spike);
    spike_y1=y1(spike_wspol)-mean(y1(spike_wspol));
    spike_y2=y2(spike_wspol)-mean(y1(spike_wspol));
    %qwer=mean(spike_y1)
    subplot(5,ls,i);
    plot(ts,spike,'b-',ts,spike_y1,'r-');
    axis([0 2.7 -200 200]);
    if i==1        
        xlabel('time [ms]');
        ylabel('signal level');
    end
    if i==2                
        legend('original signal','interpolated signal');
    end
    grid on;
    
    subplot(5,ls,i+ls);
    plot(ts,spike,'bd-',ts,spike_y1,'rd-');
    axis([0.4 1 -200 200]);
    if i==1        
        xlabel('time [ms]');
        ylabel('signal level');
    end
    if i==2        
        legend('original signal','interpolated signal');
    end
    grid on;
    
    subplot(5,ls,i+2*ls);
    fs=fft(spike,500);
    f=[0:499]/500*20000;
    plot(f,abs(fs));
    xlabel('frequecy [Hz]');
    ylabel('Fourier coefficient');
    grid on;
    %if i==1
        axis([0 10000 0 1000]);
        %end

        
%Sekcja analizy wplywu filtru z neurochipa        
    subplot(5,ls,i+3*ls);
    plot(ts,spike,'b-',ts,spike_y2,'r-');
    axis([0 1.7 -200 200]);
    if i==1        
        xlabel('time [ms]');
        ylabel('signal level');
    end
    if i==2                
        legend('original signal','interpolated signal');
    end
    grid on;
    
    subplot(5,ls,i+4*ls);
    plot(ts,spike,'bd-',ts,spike_y2,'rd-');
    axis([0.4 1 -200 200]);
    if i==1        
        xlabel('time [ms]');
        ylabel('signal level');
    end
    if i==2        
        legend('original signal','interpolated signal');
    end
    grid on;    
        
     
end
y=0;