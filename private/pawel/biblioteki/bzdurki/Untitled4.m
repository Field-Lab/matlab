clear

fid=fopen('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\dane\ID=96','r','ieee-le'); 
a=fread(fid,'int32');
l=length(a);
b=reshape(a,5,l/5);
fclose(fid);

%figure(2)
%plot(b(2,:),'bd');

%StimElectrodes0=unique(b(3,:));
%StimElectrodes=StimElectrodes0(StimElectrodes0>0);
StimElectrodes=[127 103 174 199];

NumberOfFrequencies=8;
FirstMovie=1;
NumberOfMoviesPerAmplitude=length(StimElectrodes)*NumberOfFrequencies

for i=1:length(StimElectrodes)
    StimEl=StimElectrodes(i)
    SpikesforThisEl=find(b(3,:)==StimEl);
    %MoviesForThisEl=unique(b(1,SpikesforThisEl))
    MoviesForThisEl=[[1:8]+(i-1)*8 [1:8]+(i-1)*8+32]
    
    
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
    figure(7)
    %clf
    subplot(2,2,i)
    h=plot(Xcoord,Ycoord,'rd')
    set(h,'MarkerSize',4)
    axis([0 1 0 1])
    hold on
end
    
    
