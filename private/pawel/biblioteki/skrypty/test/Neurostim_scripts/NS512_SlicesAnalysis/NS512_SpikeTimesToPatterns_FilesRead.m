clear
%'C:\pawel\nauka\analiza\slices\2010-09-14-0\data002minus009paramSninya\201
%0-09-14-0\dane\ID=6648' - dobyr neuron do testowannia multipeak fit; takze
%5823
TotalNumberOfPatterns=64;
MoviesStep=8;

neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\pawel\nauka\analiza\slices\2010-09-14-0\data002minus009paramSninya\2010-09-14-0\data002min009\data002min009.neurons');

idList = neuronFile.getIDList();
GaussParameters=double([]);

for Neuron=1:length(idList)
NeuronID=idList(Neuron)
FullName=['C:\pawel\nauka\analiza\slices\2010-09-14-0\data002minus009paramSninya\2010-09-14-0\dane\ID=' num2str(NeuronID)];;
fid=fopen(FullName,'r','ieee-le'); 
a=fread(fid,'int32');
l=length(a);
dane=reshape(a,4,l/4);
fclose(fid);

Movies=dane(1,:);
Patterns=dane(3,:);
Latencies=dane(4,:);
figure(2)
clf
p=hist(Patterns,[0:1:512]);
plot(p);

MeanSpikesPerPattern=l/4/TotalNumberOfPatterns;
MinimumNumberOfSpikesToBeInteresting=max(MeanSpikesPerPattern*2.5,190);

InterestingElectrodes=find(p>MinimumNumberOfSpikesToBeInteresting)-1; % minus1 bo pierwszy s³upek w histogramie 'p' odpowiada zeru
InterestingElectrodesFinal=InterestingElectrodes(find(InterestingElectrodes>0))

figure(1)
clf
NumberOfColumns=ceil(sqrt(length(InterestingElectrodesFinal)));
NumberOfRows=ceil(length(InterestingElectrodesFinal)/NumberOfColumns);
t=[0:1:600];
for i=1:length(InterestingElectrodesFinal)
    Electrode=InterestingElectrodesFinal(i);
    EventsForGivenElectrode=find(Patterns==Electrode);
    LatenciesForGivenElectrode=Latencies(EventsForGivenElectrode);  
    subplot(NumberOfRows,NumberOfColumns*2,(i-1)*2+1);
    p=hist(LatenciesForGivenElectrode,t);
    hist(LatenciesForGivenElectrode,t);
    h=gca;
    Yrange=get(h,'YLim');
    Ymax=Yrange(2);
    hold on;
    h=text(50,0.9*Ymax,['el. ' num2str(Electrode)]);
    set(h,'FontSize',7);
    
    CFs=NS512_FitWithMultiGauss(t,p);
    FitLine=zeros(1,length(t));
    for j=1:length(CFs)
        FitLine=FitLine+CFs{j}(t)';
        %h=text(400,(0.9-(j-1)*0.3)*Ymax,['A=' num2str(round(CFs{j}.A))]);
        %set(h,'FontSize',7);            
        h=text(320,(0.9-(j-1)*0.27)*Ymax,['\tau=' num2str(CFs{j}.tau/20,'%8.2f')]);
        set(h,'FontSize',7);
        h=text(320,(0.8-(j-1)*0.27)*Ymax,['\sigma=' num2str(CFs{j}.sigma/20,'%8.3f')]);
        set(h,'FontSize',7);
        
        DataToSave=[double(NeuronID) double( Electrode) CFs{j}.A CFs{j}.tau CFs{j}.sigma];
        GaussParameters=[GaussParameters' DataToSave']';
    end
    %gauss=cf(t);
    plot(t,FitLine,'r-');
    h=gca;
    set(h,'XLim',[0 600]);
    
    MoviesForGivenElectrode=Movies(EventsForGivenElectrode); 
    UniqueMovies=unique(MoviesForGivenElectrode);
    Eff=zeros(1,length(UniqueMovies));
    for j=1:length(UniqueMovies)
        Eff(j)=length(find(MoviesForGivenElectrode==UniqueMovies(j)));
    end
    subplot(NumberOfRows,NumberOfColumns*2,(i-1)*2+2);
    plot(1:length(UniqueMovies),Eff/108*100);
    axis([0 18 0 100]);
end
FullImageName=['C:\pawel\nauka\analiza\slices\2010-09-14-0\data002minus009paramSninya\2010-09-14-0\dane\images3\Neuron' num2str(NeuronID) '.tif'];            
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[8 4.5]);
set(h,'PaperPosition',[0 0 8 4.5]); 
print(h, '-dtiff', '-r200', FullImageName);

end
FullName=['C:\pawel\nauka\analiza\slices\2010-09-14-0\data002minus009paramSninya\2010-09-14-0\dane\images3\GaussParameters.bin'];
fid=fopen(FullName,'a','ieee-le');                                    
fwrite(fid,GaussParameters,'double');
fclose(fid);