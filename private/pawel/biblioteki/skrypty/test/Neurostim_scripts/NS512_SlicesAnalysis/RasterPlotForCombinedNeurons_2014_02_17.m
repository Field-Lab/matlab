%clear
RawDataPath='I:\analysis\slices\2013-12-12-3-PH\data005';
%paramsFilePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.params';
%neuronFilePath ='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.neurons';

paramsFile='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.params';
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\data005.neurons');
DuplicatesFilePath='D:\home\pawel\analysis\slices\2013\2013-12-12-3-PH\duplicates.txt';

NS_GlobalConstants=NS_GenerateGlobalConstants(500);
MovieFilePath='I:\analysis\slices\2013-12-12-3-PH\movie005';
%NeuronID=976;
Styles={'bo' 'ro' 'mo' 'ko'};
Offset=200;
clf
for gh=1%:length(PrimaryNeurons)
NeuronID=PrimaryNeurons(gh);
SeedEl = neuronFile.getNeuronIDElectrode(NeuronID);
    
fid=fopen(['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\dane_combined\ID=' num2str(NeuronID) 'c'],'r','ieee-le'); 
%[182 257 271 289 588 616 679 694 752 766 771 800

a=fread(fid,'int32');
l=length(a);
b=reshape(a,5,l/5); %MovieNumber,RepetitionNumber,PatternNumber,Latency,TimeRelativeToRepetitionBegin
fclose(fid);

%import direct values of spike times - necessary for EIs calculation
spikeTimes=NS512_SpikeTimesToStimulationParameters(NeuronID,paramsFile,neuronFile,MovieFilePath,DuplicatesFilePath);

StimElectrodes=[127 103 174 199];

NumberOfFrequencies=8;
FirstMovie=1;
NumberOfMoviesPerAmplitude=length(StimElectrodes)*NumberOfFrequencies;
figure(1)
clf
for i=1:length(StimElectrodes)
    SubplotID=i;%+(floor(i/3))*2;
    StimEl=StimElectrodes(i)
    SpikesforThisEl=find(b(3,:)==StimEl);
    aa=length(SpikesforThisEl)
    %MoviesForThisEl=unique(b(1,SpikesforThisEl))
    %MoviesForThisEl=[[1:8]+(i-1)*8 [1:8]+(i-1)*8+32]
    MoviesForThisEl=[[1:8]+(i-1)*8+32];    
    subplot(2,3,SubplotID)
    for movie=1:length(MoviesForThisEl)
        MovieDataFull=NS_MovieData_GlobalPath(MovieFilePath,MoviesForThisEl(movie),NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieDataFull);
        SE=MovieData(2:3:length(MovieData));
        Times=MovieData(1+(find(SE==StimEl)-1)*3);
        %subplot(2,3,SubplotID)
        x=((Times/10000)*0.9+0.05)/2;
        y=0.96-(movie-1)*0.06;
        for p=1:length(x)
            plot([x(p) x(p)],[y y-0.01],'r-');
        end
        %plot(x/2,y,'rd');
        hold on;
    end
    
    Xcoord=zeros(1,length(SpikesforThisEl));
    Ycoord=zeros(1,length(SpikesforThisEl));
    for spike=1:length(SpikesforThisEl)
        spikeID=SpikesforThisEl(spike);
        MovieNumber=b(1,spikeID);                    
            MovieID=find(MovieNumber==MoviesForThisEl)
            if length(MovieID)==0
                MovieID=-2;
            end

            AmplitudeNumber=ceil(MovieNumber/NumberOfMoviesPerAmplitude);
              
            RepetitionNumber=b(2,spikeID);
            TimeRelativeToRepetitionBegin=b(5,spikeID);
            Xcoord(spike)=TimeRelativeToRepetitionBegin/10000*0.9+0.05;
            Ycoord(spike)=0.95-((MovieID-1)*0.06+RepetitionNumber*0.0005);
    end
    
    %clf
    %subplot(2,3,SubplotID)
    h=plot(Xcoord/2,Ycoord,'bo');
    set(h,'MarkerSize',1)
    set(h,'MarkerFaceColor','b')
    %axis([0 1 0 1]) % for two amplitudes
    axis([0 0.5 0.5 1]) % for only the last amplitude
    hold on
    xlabel('Time [ms]');
    h=gca;
    set(h,'YTick',[0.94-([8:-1:1]-1)*0.06]);
    set(h,'YTickLabel',{'250' '125' '62' '31' '16' '8' '4' '2'});
    ylabel(['stimulation frequency [Hz] at electrode ' num2str(StimEl)]);
    
    %now calculate the EI but only for spikes associated with the given
    %stimulation electrode
    %SpikeTimesForElectrode=spikeTimes(SpikesforThisEl);    
    %NeuronEI=NS512_LoadNeuronEIFromRawData(RawDataPath,paramsFilePath,neuronFilePath,NeuronID,400,80,20);
    %subplot(2,3,6);
    %h=NS512_ShowEIAsCircles4(NeuronEI,500,[1:512],StimElectrodes,[],[-1000 1000],[-500 500]);
        
%end
end

subplot('position',[0.41 0.11 0.5 0.34]);
for i2=1:length(StimElectrodes)
    StimEl=StimElectrodes(i2)
    SpikesforThisEl=find(b(3,:)==StimEl);
    bb=length(SpikesforThisEl)
    if bb>50
        SpikeTimesForElectrode=spikeTimes(SpikesforThisEl);
        %NeuronEI=NS512_LoadNeuronEIFromRawData(RawDataPath,paramsFilePath,neuronFilePath,NeuronID,400,80,20);
        NeuronEI=NS512_LoadEIFromRawData(RawDataPath,SpikeTimesForElectrode,Offset+60,Offset);
        k1=max(NeuronEI')-min(NeuronEI');    
        %subplot(2,3,6);
        h=NS512_ShowEIAsEmptyCircles(NeuronEI/max(k1)*120,500,[1:512],SeedEl,StimEl,[-1000 1000],[-500 500],Styles{i2});
        hold on
    end
end


%FullImageName=['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\report_2014_02_13\figures\Neuron' num2str(NeuronID) '.tif']; % pokoj 109
FullImageName=[FigurePath 'Neuron' num2str(NeuronID) '.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]);
set(h,'PaperPosition',[0 0 16 9]); 
%print(h, '-dtiff', '-r120', FullImageName);
    
end
%trzy rzeczy:  normalizacja EI (do elektrody o najwiekszym sygnale),
%rysowac tylko EI jesli dla danej elektrody jest co najmniej 20 spikow,
%brac wartosc EI od momentu generacji impulsu stymulacyjnego do momentu
%spiku plus 1 ms?