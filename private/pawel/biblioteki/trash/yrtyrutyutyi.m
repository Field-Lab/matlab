figure(101)
for i=1:6
    subplot('position',[0.885,0.97-i*0.145,0.10,0.11]);    
    ghj=TimesForHist(i,:);
    T1=find(ghj>0);
    %T1=find(TimesForHist(i,:))~=0;
    T2=TimesForHist(i,T1);
    po=hist(T2/20,[1:1:25]/20);
    bar([1:25]/20,po/sum(po)*100,1);
    
    h=gca;
    set(h,'YLim',[0 60]);
    set(h,'YTick',[0:20:60]);
    %set(h,'YTickLabel',[0:1:3]);
    %h=ylabel('p [%]');
    set(h,'FontSize',FontSize);
    set(h,'XLim',[0.3 1.2]);
    h=ylabel('p [%]');
    
    h=text(0.8,40,['\sigma=' num2str(round(std(T2)*50)) '\mus']);
    set(h,'FontSize',FontSize);
    
    %Efektywnosc
    lk0=TimesForHist(i,:); %ilosc wystymulowanych spikow
    lk=length(find(lk0>0));
    ef=lk/IloscImpulsow(i)*100;
    h=text(0.8,20,['\epsilon=' num2str(ef,'%10.1f') '%']);
     set(h,'FontSize',FontSize);
end