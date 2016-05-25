DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc';

Electrodes=[236 244 245 412];
Electrodes=[333 334 327 328] % dla pattern 61 - axonal bundle
Electrodes=[207 189 235 218   333 334 327 328] % dla pattern 61 - pojedyncze neurony
LimityNeg=[150 200 200 300 2000 2000 2000 2000];
LimityPos=[150/2 200/2 200/2 300/2 600 600 600 600];
Limitx1=[5 2 2 0 0 0 0 0];
%Limity1=[10 7 7 0 0 0 0 0];
%Limitx1=[0 0 0 0 0 0 0 0];


Amplitudes=[1:2:25];
Pattern=61%151;
MovieStart=7;

% 229/429, 

figure(101)
clf
t=[1:400]/20;

FontSize=16;
for i=1:length(Amplitudes)
    Amplitude=Amplitudes(i);
    Movie=MovieStart+(Amplitude-1)*18;
    [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Pattern,Movie,0,0);
    
    for j=1:length(Electrodes)
        data=reshape(DataTraces(:,Electrodes(j),:),50,400);
        %subplot(length(Electrodes),length(Amplitudes),length(Amplitudes)*(j-1)+i);
        subplot('position',[(i-1)*0.07+0.07 (9-j)*0.11  0.06 0.08]);   

        h1=plot(t,data'/0.27);
        if j<5
            set(h1,'Color','r')
        else
            set(h1,'Color','b');
        end
        axis([Limitx1(j) Limitx1(j)+5 -LimityNeg(j) LimityPos(j) ])
        grid on
        h2=gca
        set(h2,'XTick',[Limitx1(j):1:Limitx1(j)+5]);        
        set(h2,'FontSize',FontSize);
        
        if j==length(Electrodes) && i==1
            h3=xlabel('Time [ms]');
            set(h3,'FontSize',FontSize);
            h3=ylabel('Signal [\muV]');
            set(h3,'FontSize',FontSize);
            %ylabel('Signal [\muV]');
        end
        if i==1
            %h3=ylabel('Signal [\muV]');
            %set(h3,'FontSize',FontSize);
        else
            set(h2,'XTickLabel','');
            set(h2,'YTickLabel','');
        end
    end
end

FullName=[FigurePath '\Figure3_E.tif']; 
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 10]);
set(h,'PaperPosition',[0 0 16 10]); 
print(h, '-dtiff', '-r100', FullName);
