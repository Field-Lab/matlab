NS_GlobalConstants=NS_GenerateGlobalConstants(500);
% 1) Find primary electrode for each cell
NeuronIDs=[6 110 139 320 349 367 636 726 873 889 1026 1190 1283 1428 1445 1518 1732 1804 1880 1908 2017 2044 2089 2269 2299 2483 2677 2722 2900 3037 3218 3245 3425 3605 3622 3917 4055 4158 4310 4401 4578 4609 4865 4895 4907 5102 5162 5255 5389 5602 5633 5798 5884 5977 6064 6110 6263 6320 6425 6473 6503 7011 7025 7069 7249 7278 7324 7427 7461 7638];

%NeuronIDs=[6 110 139 320 349 367 636 726 873 889 1026 1190 1283 1428];

RawDataPath='E:\2012-09-27-4\data000';
VisionOutputPath='E:\analysis\2012-09-27-4\data000';
%VisionOutputPath='D:\Home\Pawel\analysis\2012_09\2012-09-27-4\data000';
PrimaryElectrodes=NeuronIDs;
EIs=zeros(length(NeuronIDs),60);

for i=1:length(NeuronIDs)
    ID=NeuronIDs(i);
    [MinValues,MaxValues,EI]=NS512_FindEIamplitudes(RawDataPath,VisionOutputPath,'000',ID);    
    Amplitudes=MaxValues-MinValues;
    [A,I]=max(Amplitudes);
    PrimaryElectrodes(i)=I;
    MaxY=MaxValues(PrimaryElectrodes(i));
    MinY=MinValues(PrimaryElectrodes(i));
    HighLimit(i)=MaxY+(MaxY-MinY)*0.1;
    LowLimit(i)=MinY-(MaxY-MinY)*0.1;  
    
    EIs(i,:)=reshape(EI(1,I,:),1,60);
end

%2) Load information about the stimulating electrode and amplitude for each
%neuron
ThresholdFilePath1=['E:\analysis\2012-09-27-4\stim_scan\thresholds_1']
fid1=fopen(ThresholdFilePath1,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
NeuronInformation=reshape(b,length(b)/4,4);
Neurons1=NeuronInformation(:,1);
Electrodes1=NeuronInformation(:,3);
Amplitudes1=NeuronInformation(:,4);

ThresholdFilePath2=['E:\analysis\2012-09-27-4\stim_scan\thresholds_2']
fid1=fopen(ThresholdFilePath2,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
NeuronInformation=reshape(b,length(b)/4,4);
Neurons2=NeuronInformation(:,1);
Electrodes2=NeuronInformation(:,3);
Amplitudes2=NeuronInformation(:,4);

Neurons=[Neurons1' Neurons2']';
Electrodes=[Electrodes1' Electrodes2']';
Amplitudes=[Amplitudes1' Amplitudes2']';
%trzeba napisac funkcje przeksztalcajaca amplitude na numer movie

%3) Amplitudes -> movie numbers
%StimulationDataPath='D:\Home\Pawel\analysis\2012_09\2012-09-27-4\scan_new';
StimulationDataPath='E:\analysis\2012-09-27-4\scan_new';

%FigurePath='D:\Home\Pawel\analysis\2012_09\2012-09-27-4\figures';
%FigurePath='E:\analysis\2012-09-27-4\figures';
FigurePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\figures_try';

Movies=[1:2:63];

for i=1:length(Movies)
    MovieNumber=Movies(i);
    %MoviesAmplitudes(i)=NS_AmplitudesForPattern_512_1el('D:\Home\Pawel\analysis\2012_09\2012-09-27-4\scan_new',2,2,MovieNumber,NS_GlobalConstants);
    MoviesAmplitudes(i)=NS_AmplitudesForPattern_512_1el('E:\analysis\2012-09-27-4\scan_new',2,2,MovieNumber,NS_GlobalConstants);
end
    
for i=1:length(Amplitudes)
    [y,MovieIndex]=min(abs(MoviesAmplitudes-Amplitudes(i)));
    ThresholdMovies(i)=Movies(MovieIndex);
end

NumberOfColumns=10;
NumberOfRows=7;
Patterns=Electrodes;
FontSize=10;
Crosstalk=zeros(length(Neurons),length(Neurons));
for i=1:length(Patterns)
    Pattern=Patterns(i);
    StimulatedNeuron=Neurons(i)
    %for j=1:length(Movies)
        Movie=ThresholdMovies(i);
        [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(StimulationDataPath,StimulationDataPath,0,Pattern,Movie,0,0);   %read data                 
        SDT0=size(DataTraces0);
        DataTraces=DataTraces0(:,PrimaryElectrodes,8:SDT0(3));
        SDT=size(DataTraces);
        figure(1);
        clf;
        for k=1:length(PrimaryElectrodes)
            RecordedNeuron=NeuronIDs(k)
            s=reshape(DataTraces(:,k,:),SDT(1),SDT(3));
            RowIndex=ceil(k/NumberOfColumns);
            ColumnIndex=k-(RowIndex-1)*NumberOfColumns;
            subplot('Position',[0.02+(ColumnIndex-1)*0.98/NumberOfColumns 0.99-RowIndex*0.98/NumberOfRows 0.98/NumberOfColumns-0.02 0.98/NumberOfRows-0.02]);
                                    
            h6=plot(EIs(k,:));
            set(h6,'Color','b');
            set(h6,'LineWidth',2);
            set(h6,'LineStyle','--'); 
            hold on;
            offset=EIs(k,1);
            minimum=min(EIs(k,:));
            CrosstalkThreshold=offset+(minimum-offset)*0.7;
            CrosstalkSpikes=find(min(s')-CrosstalkThreshold<0); %at which stimulation pulses repetitions, the given (non-stimulated) primaty electrode shows signals that may be spikes
            Crosstalk(i,k)=length(CrosstalkSpikes);
            h7=text(100,HighLimit(k)-(HighLimit(k)-LowLimit(k))*0.1,num2str(length(CrosstalkSpikes)));
            set(h7,'FontSize',FontSize);
                       
            h5=plot(s');
            for l=1:SDT(1)
                if find(CrosstalkSpikes==l) & length(CrosstalkSpikes)>5
                    set(h5(l),'Color','r');
                else
                    set(h5(l),'Color','k');
                end
            end                
            %set(h5,'Color','k');
            h1=gca;
            set(h1,'YTickLabel',[]);
            set(h1,'FontSize',FontSize);
            set(h1,'XTickLabel',[]);            
            set(h1,'XLim',[0 120]);
            set(h1,'YLim',[LowLimit(k) HighLimit(k)]);
            h2=text(5,HighLimit(k)-(HighLimit(k)-LowLimit(k))*0.1,[num2str(PrimaryElectrodes(k))]);
            set(h2,'FontSize',FontSize);
            if StimulatedNeuron==RecordedNeuron               
                set(h5,'Color','g');
            end            
            %hold on;
            subplot('Position',[0.02+(ColumnIndex-1+0.5)*0.98/NumberOfColumns 0.99-(RowIndex+0.7)*0.98/NumberOfRows (0.98/NumberOfColumns)*0.4-0.02 (0.98/NumberOfRows)*0.2-0.02]);
        end
        FigureName=[FigurePath '\p' num2str(Pattern) 'm' num2str(Movie)];            
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]);  
        print(h, '-dtiff', '-r120', FigureName);
    %end
end