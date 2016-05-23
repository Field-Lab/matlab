function InterestingElectrodes=NS512_SpikesAnalysis_FindElectrodesWithTimeLockedResponses(spfilepath,Pattern,Movies,TimeWindow);
% Najpierw histogramy spikow dla wszystkich movies s¹ sumowane (2014/12/02)

SpikesPerElectrode=zeros(length(Movies),512);
AllDelaysHistograms=zeros(512,600);
for m=1:length(Movies)
    Movie=Movies(m);
    [spfilepath 'sp_p' num2str(Pattern) 'm' num2str(Movie)]
    fid=fopen([spfilepath 'sp_p' num2str(Pattern) 'm' num2str(Movie)],'r','ieee-le');
    a=fread(fid,'int32');
    b=reshape(a,length(a)/3,3);
    fclose(fid);
    
    for electrode=1:512
        %size(AllDelaysHistograms(electrode,:))
        %size(hist(b(find(b(:,1)==electrode),3),[1:600]))
        h1=hist(b(find(b(:,1)==electrode),3),[1:600]);
        h2=reshape(h1,1,600);
        AllDelaysHistograms(electrode,:)=AllDelaysHistograms(electrode,:)+h2;
        %AllDelaysHistograms(electrode,:)=AllDelaysHistograms(electrode,:)+hist(b(find(b(:,1)==electrode),3),[1:600]);
    end
        
    SpikesPerElectrode(m,:)=hist(b(:,1),[1:512]);
end

for i=1:512
    TimeLocking(i)=NS512_AreTheSpikesTimeLocked(AllDelaysHistograms(i,:),TimeWindow); % na tych elektrodach dzieje sie *cos" ciekawego
end

InterestingElectrodes=find(TimeLocking==1);