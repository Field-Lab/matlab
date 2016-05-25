NS_GlobalConstants=NS_GenerateGlobalConstants(500);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);

GaussesFilesPath='D:\Home\Pawel\analysis\slices\512stim_paper\SpikesAnalysis\2010-09-14-0\GaussesFiles2\';

for i=1:512
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);        
end

NumberOfAmplitudes=17;
NumberOfStimulatedElectrodes=zeros(NumberOfAmplitudes,512);
AllPatterns=[]
figure(11)
clf
 for MovieSeqeunce=1:NumberOfMovieSequences
        Movie=(Amplitude-1)*NumberOfMovieSequences2+MovieSeqeunce;
        ChunkData=NS_MovieData_GlobalPath(MovieFile,Movie,NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);  
        
        PatternsForMovie=unique(MovieData(2:3:length(MovieData)));
        AllPatterns=[AllPatterns' PatternsForMovie']';
 end
 
 
 
        for p=1:length(AllPatterns)
            Pattern=AllPatterns(p);            
            %n=StimulatedElectrodesMouse(Amplitude,Pattern,:);
            %ElectrodesWithSpikes=find(n>25);
            fid=fopen([GaussesFilesPath 'AllGausses_p' num2str(Pattern)],'r');
            if fid~=-1
            a=fread(fid,'double');
            fclose(fid);    
            b=reshape(a,length(a)/5,5); 
            
            ElectrodesWithConnection0=b(:,2);
            ElectrodesWithLateSpikes=find(b(:,4)>0);
            ElectrodesWithConnection=ElectrodesWithConnection0(ElectrodesWithLateSpikes);
            
            clear X;
            clear Y;
            clear X2;
            clear Y2;
            for i=1:length(ElectrodesWithConnection)
               X(i)=electrodeMap.getXPosition(ElectrodesWithConnection(i));
               Y(i)=electrodeMap.getYPosition(ElectrodesWithConnection(i));
            end
            
            PatternX=electrodeMap.getXPosition(Pattern);
            PatternY=electrodeMap.getYPosition(Pattern);
            
            SubplotColumn=(PatternX+945)/120;
            SubplotRow=(PatternY+450)/120;
            
            subplot('position',[SubplotColumn*0.06+0.01 SubplotRow*0.06+0.02 0.054 0.048]);   
            
            if ElectrodesWithConnection
                h1=plot(X,Y,'ro');
                set(h1,'MarkerFaceColor','r');
                set(h1,'MarkerSize',4)
            end
            hold on
            
            h1=plot(PatternX,PatternY,'gd');    
            set(h1,'MarkerFaceColor','g');
            axis([-1000 1000 -500 500]);
            h=gca;
            set(h,'XTickLabel','');
            set(h,'YTickLabel',''); 
            set(h,'Box','on');
            %h=text(450,-360,num2str(Pattern));
            %set(h,'FontSize',8); 
            end
                                                                                                        
        end       
%    end        
%end