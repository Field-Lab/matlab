%clear

NS_GlobalConstants=NS_GenerateGlobalConstants(500);
MovieFilePath='I:\analysis\slices\2013-12-12-3-PH\movie005';
%NeuronID=976;

for gh=1:3
NeuronID=PrimaryNeurons(gh);
    
fid=fopen(['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\dane_combined\ID=' num2str(NeuronID) 'c'],'r','ieee-le'); 
%[182 257 271 289 588 616 679 694 752 766 771 800

a=fread(fid,'int32');
l=length(a);
b=reshape(a,5,l/5);
fclose(fid);

%figure(2)
%plot(b(2,:),'bd');

%StimElectrodes0=unique(b(3,:));
%StimElectrodes=StimElectrodes0(StimElectrodes0>0);
StimElectrodes=[127 103 174 199];
%StimElectrodes=[127 103 174 199];

NumberOfFrequencies=8;
FirstMovie=1;
NumberOfMoviesPerAmplitude=length(StimElectrodes)*NumberOfFrequencies
figure(1)
clf
for i=1:length(StimElectrodes)
    StimEl=StimElectrodes(i)
    SpikesforThisEl=find(b(3,:)==StimEl);
    %MoviesForThisEl=unique(b(1,SpikesforThisEl))
    %MoviesForThisEl=[[1:8]+(i-1)*8 [1:8]+(i-1)*8+32]
    MoviesForThisEl=[[1:8]+(i-1)*8+32]    
    for movie=1:length(MoviesForThisEl)
        MovieDataFull=NS_MovieData_GlobalPath(MovieFilePath,MoviesForThisEl(movie),NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieDataFull);
        SE=MovieData(2:3:length(MovieData));
        Times=MovieData(1+(find(SE==StimEl)-1)*3);
        subplot(2,2,i)
        x=((Times/10000)*0.9+0.05)/2;
        y=0.96-(movie-1)*0.06;
        for p=1:length(x)
            plot([x(p) x(p)],[y y-0.01],'r-');
        end
        %plot(x/2,y,'rd');
        hold on;
    end
        
    for spike=1:length(SpikesforThisEl)
        spikeID=SpikesforThisEl(spike);
        MovieNumber=b(1,spikeID);                    
            MovieID=find(MovieNumber==MoviesForThisEl);
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
    subplot(2,2,i)
    h=plot(Xcoord/2,Ycoord,'bo')
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
end

FullImageName=['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\report_2014_02_13\figures\Neuron' num2str(NeuronID) '.tif']; % pokoj 109
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[8 8]);
set(h,'PaperPosition',[0 0 8 8]); 
print(h, '-dtiff', '-r240', FullImageName);
    
end