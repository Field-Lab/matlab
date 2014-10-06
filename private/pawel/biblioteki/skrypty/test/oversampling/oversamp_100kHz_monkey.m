%10:13 dla cell body
%1,2,3,5 dla axonal
system=2; %windows
figura1=5;
figura2=6;
figura3=105;
figura4=106;
figura5=205;
figura7=305;

fft_marg=10; % wartosc marginesu dla funkcji fft_blackman 
N_fft=400;
fp=20000;

ktore_spiki=[1:16];
%ktore_spiki=[10:13];
%ktore_spiki=[1 4 5 15];
[p,filenames,spikes,detect_param]=spiki_malpa(1,ktore_spiki,system);
%[p,filenames,spikes,detect_param]=spiki_gunpg(1,ktore_spiki,system);

dl=400;
[y1,y2]=inveyefilter_hayes(dl,58,40,3300,20000);

%
[y1a,y2a]=inveyefilter_hayes(dl,18,20,8300,20000);

margins=[15 50];
t=[0:margins(1)+margins(2)]/20;
a1=4;
a2=4; % ile wykresow - wpsolrzedne dla subplot
%*  *  *  *  * 1 - rysowanie 16 spikow type cell body, duza skala czasowa
%figure(5);
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
            figure(figura1);
	        subplot(a1,a2,positions(i));
            if positions(i)==1
                lgnd=1;
            else
                lgnd=0;
            end
            n=oversamp_100kHz_plot(signal,10,lgnd,0);            	        		    
		%h=gca;
		%set(h,'XLim',[0 3]);
		%xlabel('time [ms]');         
        end
    end    
end

[p,filenames,spikes,detect_param]=spiki_malpa(2,ktore_spiki,system);
ktore_spiki=[1:16];
%ktore_piki=[2 3 5 12];
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
            figure(figura2);
	        subplot(a1,a2,positions(i));
            if positions(i)==1
                lgnd=1;
            else
                lgnd=0;
            end
            n=oversamp_100kHz_plot(signal,10,lgnd,0);            	        		    
		%h=gca;
		%set(h,'XLim',[0 3]);
		%xlabel('time [ms]');         
        end
    end    
end



%[p,filenames,spikes,detect_param]=spiki_malpa(2,ktore_spiki,system);

%dl=200;
%[y1,y2]=inveyefilter_hayes(dl,58,40,3300,20000);
