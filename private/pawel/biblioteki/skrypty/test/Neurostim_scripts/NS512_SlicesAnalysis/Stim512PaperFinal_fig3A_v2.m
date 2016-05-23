NS_GlobalConstants=NS_GenerateGlobalConstants(500);
[NumberOfMovies,NumberOfPatternsPerMovie,AllPatterns,Patterns]=NS512_MoviePatterns(MovieFile,NS_GlobalConstants);

StimulatedElectrodes=StimulatedElectrodesMouse;
s1=size(StimulatedElectrodes);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
SPFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\sp_files\';

for i=1:512
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);        
end

MapaKolorow=jet;
clear qw;
PatternsToAnalyze=unique(AllPatterns)';
KtoryPattern=17;



for m=1:0%NumberOfAmplitudes %%% czasem od 2!!!
            Movie=Movies(m);
            fid=fopen([SPFilesPath 'sp_p' num2str(Pattern) 'm' num2str(Movie)],'r','ieee-le');
            a=fread(fid,'int32');
            b=reshape(a,length(a)/3,3);
            fclose(fid);    
        
            SpikesIDs=find(b(:,1)==Electrode);

            Delays=b(SpikesIDs,3);
            AllDelays=[AllDelays' Delays']';           
        end


SpikesOverall=zeros(1,512);
for i=KtoryPattern
    Pattern=PatternsToAnalyze(i)
    %clf
    %StimulatedChannel=PatternsToAnalyze(i);
    %ChannelData=StimulatedElectrodes(:,StimulatedChannel,:);
    %qw(StimulatedChannel)=max(max(max(ChannelData)));
    
    %ChannelDataNew=reshape(ChannelData,s1(1),s1(3));
    ChannelsToExclude=electrodeMap.getAdjacentsTo(Pattern,1)';
    %ChannelDataNew(:,ChannelsToExclude)=0; %ignore spikes on the stimulated electrode
    
    %qw2(StimulatedChannel)=max(max(max(ChannelData)));
    
    for Amplitude=2:25
        Movie=7+(Amplitude-1)*18;
        fid=fopen([SPFilesPath 'sp_p' num2str(Pattern) 'm' num2str(Movie)],'r','ieee-le');
        a=fread(fid,'int32');
        b=reshape(a,length(a)/3,3);
        fclose(fid);  
        
        Spikes=zeros(1,512);
        for el=1:512
            l=find(b(:,1)==el);
            if l
                Spikes(el)=length(l);
            end
        end
        Spikes(ChannelsToExclude)=0;
        SpikesOverall=max(SpikesOverall,Spikes);
        %Amplitude
        SubplotRow=ceil((Amplitude-1)/6)
        SubplotColumn=Amplitude-1-(SubplotRow*6)+6-1;
        subplot('position',[SubplotColumn*0.12+0.01 (4-SubplotRow)*0.1+0.56 0.1125 0.09]);        
        
        %ChannelDataForAmplitude=ChannelDataNew(Amplitude,:);
        ConvertedData=ChannelDataToColors(Spikes,X,Y);
        image(ConvertedData');
        colormap(jet(50))
        hold on
        h10=plot(X(Pattern)/30+33,17-round((Y(Pattern)+490)/60),'go');
        set(h10,'MarkerSize',8)
        %axis([-1000 1000 -500 500])
        h=gca;        
        %set(h,'Color',[0 0 0]);
        set(h,'XTickLabel','');
        set(h,'YTickLabel',''); 
        %et(h,'Box','on');
    end    
end

