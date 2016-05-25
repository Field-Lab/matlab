clear

% do ustalenia: a) najmniejsza wartoœæ czêstotliwoœci dla prac: Jankowska1972,
% Jankowska1973; Stoney1968, Tehovnik2009; b) powierzchnia elektrody dla Tehovnik2003, Tehovnik2009,
% Houweling2008, Murphey2007; c)
% dodac wyniki dla stymulacji nerwu s³uchowego 

%inne papier: Breckenridge1995 - nie podaj¹ czasu trwania impulsu;
%Gustafsson1976 - nie podaj¹ czêstotliwoœci impulsów; 

%lit(1)=struct('name','Butovas et al, 2003', 'number', '1', 'current', '8', 'charge', '800', 'area', '14130', 'f', '1','detect','extra','region','cortex');
lit(1)=struct('name','Gunning et al, 2013', 'number', '1', 'current', '1.1', 'charge', '55', 'area', '78', 'f', '2','detect','extra','region','ctx'); % area do korekty - nspisany mail do Debbie
%lit(2)=struct('name','Histed et al, 2009', 'number', '5', 'current', '4', 'charge', '800', 'area', '12.56', 'area_min', '7.07', 'area_max', '19.63', 'f', '250','detect','calcium','region','ctx'); % wyniki dla pipety szklanej, elektroda metalowa duzo wieksza ale prad nadal tylko 8 uA
lit(2)=struct('name','Histed et al, 2009', 'number', '2', 'current', '4', 'charge', '800', 'area', '12.56', 'f', '250','detect','calcium','region','ctx'); % wyniki dla pipety szklanej, elektroda metalowa duzo wieksza ale prad nadal tylko 8 uA
lit(3)=struct('name','Houweling & Brecht, 2008', 'number', '3', 'current', '2', 'charge', '600', 'area', '141300000', 'f', '200','detect','behav','region','ctx'); % œrednica pipety szklanej nie jest podana ani w artykule, ani w suplementary
lit(4)=struct('name','Murphey & Maunsell, 2007', 'number', '4', 'current', '3', 'charge', '600', 'area', '141300000', 'f', '200','detect','behav','region','ctx'); % wielkosc elektrody niepodana, tylko impedancja 0.2-1.5 MOhm
lit(5)=struct('name','Nowak & Bullier, 1998', 'number', '5', 'current', '2.5', 'charge', '250', 'area', '700', 'f', '0.5','detect','intra','region','ctx'); % niska czestotliwosc
lit(6)=struct('name','Stoney et al, 1968', 'number', '6', 'current', '1', 'charge', '200', 'area', '122.65', 'f', '-1','detect','extra','region','ctx'); % troche zamieszania z typem elektrody, uznajemy ze wire electrode 12.5 um
lit(7)=struct('name','Tehovnik et al, 2003', 'number', '7', 'current', '3', 'charge', '600', 'area', '141300000000', 'f', '200','detect','nie wiem','region','ctx'); % brak info o roymiarye elektrod
lit(8)=struct('name','Tehovnik & Slocum, 2009', 'number', '8', 'current', '6', 'charge', '1200', 'area', '1413000000000', 'f', '-1','detect','behav','region','ctx'); % brak info o roymiarye elektrod
lit(9)=struct('name','Wagenaar et al, 2004', 'number', '9', 'current', '1', 'charge', '160', 'area', '706.5', 'f', '1','detect','extra','region','ctx'); % drugi wynik: 50 us, 2 uA
%lit(10)=struct('name','Swadlow, 1992', 'number', '10', 'current', '0.2', 'charge', '60', 'area', '1413000000000', 'f', '1','detect','extra','region','ctx'); % drugi wynik: 50 us, 2 uA

lit_sp(1)=struct('name','Jankowska & Roberts, 1972', 'number', '11', 'current', '0.1', 'charge', '10', 'area', '3.14', 'f', '300','detect','extra','region','spc'); % nie wiem co z czestotliwoscia
lit_sp(2)=struct('name','Jankowska & Smith, 1973', 'number', '12', 'current', '0.6', 'charge', '30', 'area', '3.14', 'f', '300','detect','extra','region','spc'); % nie wiem czy nie wywalic?
lit_sp(3)=struct('name','Gustafsson & Jankowska, 1976', 'number', '13', 'current', '0.15', 'charge', '30', 'area', '3.14', 'f', '-1','detect','extra','region','spc'); % OK

lit_ret(1)=struct('name','Sekirnjak et al, 2006', 'number', '14', 'current', '0.6', 'charge', '60', 'area', '28.26', 'f', '1','detect','extra','region','ret'); % OK
lit_ret(2)=struct('name','Jepson et al, 2013', 'number', '15', 'current', '0.3', 'charge', '15', 'area', '176', 'f', '1','detect','extra','region','ret'); % OK
lit_ret(3)=struct('name','Jensen et al, 2003', 'number', '16', 'current', '0.26', 'charge', '26', 'area', '34.44', 'f', '4','detect','not sure','region','ret');
lit_ret(4)=struct('name','Grumet et al, 2003', 'number', '17', 'current', '0.06', 'charge', '24', 'area', '78.5', 'f', '-1','detect','not sure','region','ret');

data(1)=struct('name','Hottowy et al, 2015', 'number', 'this work', 'current', '0.3', 'charge', '30', 'area', '19.62', 'f', '2','detect','extra','region','ctx');

ll=length(lit);
for i=1:length(lit_sp)
    lit(ll+i)=lit_sp(i);
end

ll=length(lit);
for i=1:length(lit_ret)
    lit(ll+i)=lit_ret(i);
end

ll=length(lit);
for i=1:length(data)
    lit(ll+i)=data(i)
end

MS=12;
FontSize1=12

figure(1)
clf
lit2=lit;
clear lit;
lit(1)=struct('name','Murphey & Maunsell, 2007', 'number', '7', 'current', '300', 'charge', '60000', 'area', '14130', 'f', '200','detect','behav','region','ctx');
lit(2)=struct('name','Murphey & Maunsell, 2007', 'number', '7', 'current', '300', 'charge', '60000', 'area', '14130', 'f', '200','detect','behav','region','spc');
lit(3)=struct('name','Murphey & Maunsell, 2007', 'number', '7', 'current', '300', 'charge', '60000', 'area', '14130', 'f', '200','detect','behav','region','ret');
lit(4:3+length(lit2))=lit2;
for i=1:length(lit)
    name=lit(i).name
    current=lit(i).current    
    charge=lit(i).charge
    f=lit(i).f
    detect=lit(i).detect
    area=lit(i).area
    number=lit(i).number;
    %number=num2str(i-3)
    region=lit(i).region;
    detect=lit(i).detect;
    
    label=['[' number ']'];
    
    switch region
        case 'ctx'
            symc='k';
        case 'spc'
            symc='b'
        case 'ret'
            symc='r'
        otherwise
            symc='g';
    end
    if i==length(lit)
        symc='k';
        label='this work';
    end       
    %figure(1)
    subplot(1,2,1)
    hold on
    h=plot(str2num(current),str2num(f),'bd');
     set(h,'Markersize',MS');
    set(h,'MarkerEdgeColor',symc);
    set(h,'MarkerFaceColor',symc);
    h1=text(str2num(current)*1.05,str2num(f),label)
    set(h1,'FontSize',FontSize1);
    xlabel('current [uA]');
    ylabel('frequency [Hz]')
        
    subplot(1,2,2)
    hold on
    h=plot(str2num(charge),str2num(charge)/str2num(area)*1e-8,'bd');
    set(h,'Markersize',MS');
    set(h,'MarkerEdgeColor',symc);
    set(h,'MarkerFaceColor',symc);

    h1=text(str2num(charge)*0.9,str2num(charge)/str2num(area)*1.33e-8,label)
    set(h1,'FontSize',FontSize1);
    xlabel('charge [pC]');
    ylabel('charge density [pC/cm^2]')
    
    %figure(2)
    
    %area2=sqrt(area);
    %subplot(1,3,3)
    %hold on
    %h=plot(str2num(charge),str2num(charge)/sqrt(str2num(area))*1e-8,'bd');
    %h=plot(str2num(charge),str2num(area),'bd');
    %set(h,'Markersize',MS');
    %set(h,'MarkerEdgeColor',symc);
    %set(h,'MarkerFaceColor',symc);

    %h1=text(str2num(charge)*0.9,str2num(charge)/sqrt(str2num(area))*1.33e-8,label)
    %set(h1,'FontSize',FontSize1);
    %xlabel('charge [pC]');
    %ylabel('charge density [pC/cm^2]')
end

subplot(1,2,1)
h=gca
set(h,'YScale','log')
set(h,'XScale','log')
set(h,'FontSize',FontSize1);
grid on
axis([0.05 8 0.1 1000])
subplot(1,2,2)
h=gca
set(h,'YScale','log')
set(h,'XScale','log')
set(h,'FontSize',FontSize1);

grid on
axis([5 2000 1e-10 1e-6])
legend('cortex','spinal cord','retina','Location','NorthWest')
%subplot(1,3,3)
%h=gca
%set(h,'YScale','log')
%set(h,'XScale','log')
%set(h,'FontSize',FontSize1);

%grid on
%axis([10 300 1 1000])

break
subplot(2,2,3)
h=gca
set(h,'YScale','log')
set(h,'XScale','log')
grid on
subplot(2,2,4)
h=gca
set(h,'YScale','log')
set(h,'XScale','log')
grid on
%h=gca;
%set(h,'YScale','log')
%axis([0 10 0.1 300])
%grid on