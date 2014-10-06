%addpath /home/pawel/pliki/matlab/tulboksy/biblioteki/funkcje/test/stymulacja/przekladki;
%cd H:\pliki\nauka\stymulacja\chip\symulacje;
%2005_01_24: symulacje w spectre, 2005_01_27: w hspice

%cd /home/pawel/pliki/stymulacja/chip/symulacje/2005_01_24;
currents={'1m' '250u' '60u' '15u' '4u' '1u' '250n' '60n'};
res={'1.5k' '6k' '24k' '96k' '380k' '1.5M' '6M' '24M'};

%cd /home/pawel/pliki/stymulacja/chip/symulacje/2005_04_01; %symulacje w spectre
cd H:\pliki\nauka\stymulacja\chip\symulacje\2005_04_01;
res={'1.45k' '5.8k' '23k' '90k' '360k' '1.44M' '5.8M' '23M'};
res={'720' '2.9k' '11.5k' '45k' '180k' '720k' '2.9M' '11.5M'};
%res={'shorted' 'shorted' 'shorted' 'shorted' 'shorted' 'shorted' 'shorted' 'shorted'};

figure(21);
clf;
for i=1:8
    subplot(2,4,i);
    cr=currents(i)
    r=res(i);
    if i==5
	    osie=1;
    else
	    osie=0;
    end
    y=rysuj_bledy3(cr{1},r{1},7,osie);
    %y=rysuj_bledy3(cr{1},'shorted',7,osie);
    %min(y)
    %max(y)
    axis([-127 127 -1 1]);
    h=gca;
    set(h,'FontSize',16);

        
    %pos/neg
end

figure(22);
clf;
for i=1:8
    subplot(2,4,i);
    cr=currents(i)
    r=res(i);
    if i==5
	    osie=1;
    else
	    osie=0;
    end
    y=symetria1(cr{1},r{1},7,osie);
    %y=symetria1(cr{1},'shorted',7,osie);
    %min(y)
    %max(y)
    axis([0 127 0.98 1.02]);
    h=gca;
    set(h,'FontSize',16);
    %axis([32 127 0.99 1.01]);
        
    %pos/neg
end
%print -dtiff -r300 figura

%y=Rout1('15u','shorted','90k',7,1);
figure(23);
clf;
for i=1:8
    subplot(2,4,i);
    cr=currents(i)
    r=res(i)
    if i==5
	    osie=1;
    else
	    osie=0;
    end
    y=Rout1(cr{1},'shorted',r{1},1,7,osie);
    %min(y)
    %max(y)
    axis([0 127 1e1 1e5]);
    h=gca;
    set(h,'FontSize',16);
    %axis([32 127 0.99 1.01]);
        
    %pos/neg
end
