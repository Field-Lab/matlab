function y=oversamp_100kHz_report(animal,spike_type,figura,ktore_spiki,style);
%animal: 0 - malpa, 1 - guinepig
%spike_type: 0 - cell body, 1 - other
%figura - numer rysunku
%ktore_spiki - wiadomo
%style: 0 - duza skala czasowa, 1 - details (dla kazdego spikea dwa
%wykresy)

system=1; %windows
fp=20000;

if(animal==0)
    [p,filenames,spikes,detect_param]=spiki_malpa(spike_type,ktore_spiki,system);
else
    [p,filenames,spikes,detect_param]=spiki_gunpg(spike_type,ktore_spiki,system);
end

for j=1:length(filenames)
    if find(spikes(:,1)==j)'
        s=importdata(filenames{j})';
        s=s-mean(s);
    
        positions=find(spikes(:,1)==j)' %numery dobrych spikow w tym pliku    
        %coordinates(positions)=y(1,spikes(positions,2)); 

        for i=1:length(positions)
            strt=max(spikes(positions(i),2)-2000,1);
	        stp=min(spikes(positions(i),2)+2000,length(s));
            signal=s(1,strt:stp);
            figure(figura);
            
            if(style==0)
                a1=4;
                a2=4;
                subplot(a1,a2,positions(i));
                if positions(i)==1
                    lgnd=1;
                else
                    lgnd=0;
                end
                n=oversamp_100kHz_plot(signal,10,lgnd,0); 
            else
                subplot(2,length(ktore_spiki),positions(i));        
                if positions(i)==1
                    lgnd=1;
                else
                    lgnd=0;
                end
                n=oversamp_100kHz_plot(signal,10,lgnd,1);   
                h=gca;
                set(h,'XLim',[99 103]);
            
                subplot(2,length(ktore_spiki),positions(i)+length(ktore_spiki));            
                n=oversamp_100kHz_plot(signal,10,lgnd,0)    
                h=gca;
                set(h,'XLim',[n-0.05 n+0.05]);
            end
		%h=gca;
		%set(h,'XLim',[0 3]);
		%xlabel('time [ms]');         
        end
    end    
end
y=1;
