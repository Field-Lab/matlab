%cd H:\pliki\nauka\stymulacja\chip\symulacje;
cd /home/pawel/pliki/stymulacja/chip/symulacje/2005_01_24;

currents={'1m' '250u' '60u' '15u' '4u' '1u' '250n' '60n'};
res={'1.5k' '6k' '24k' '96k' '380k' '1.5M' '6M' '24M'};

figure(1);
clf;
for i=1:8
    subplot(2,4,i);
    cr=currents(i)
    r=res(i);
    y=rysuj_bledy2(cr{1},r{1},7);
    axis([0 1 -2 2]);
    
    filename=[cr{1} '_neg_' r{1} '.out'];
    a=importdata(filename);
    x=a(:,1);
    y=a(:,2);
    min(y)
    neg=y(round(length(y)/2));
    
    filename=[cr{1} '_pos_' r{1} '.out'];
    a=importdata(filename);
    x=a(:,1);
    y=a(:,2);
    max(y)
    pos=y(round(length(y)/2));
    
    %pos/neg
end



%print -dtiff -r300 figura