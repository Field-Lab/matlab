function y=spikes_oversamp_stat1(s,detect_param,filter_param,margins,figures,RC_filter);

s(1,1:margins(1,1))=0;
s(1,length(s):-1:length(s)-margins(1,2))=0;
y=blad_vs_energy2(s,detect_param,filter_param,margins);

%1. Energia bledu vs energia spika:
figure(figures(1,1));
plot(y(2,:),y(3,:),'b.');

%2. Plotowanie 300 spikow:
figure(figures(1,2));
for i=1:min(length(y),300)
    subplot(15,20,i);
    plot(s(y(1,i)-margins(1,1):y(1,i)+margins(1,2)));
    axis([0 sum(margins) -200 200]);
end

%3. Wartosci bledu dla kolejnych spikow - przydatne do znalezienia numeru
%spiku o np. najwiekszej wartosci bledu:
figure(figures(1,3));
plot(y(3,:));

%4. Blizsza analiza kilku wybranych spikow:
spikes=[35 66 64 197 207];
y=oversamp_report_onechannel(s,y(1,:),spikes,filter_param,figures(1,4),2,margins,RC_filter);
