Movie=223;
NumberOfElectrodes=4;
NumberOfRepetitions=50;
Length=10000;

StimElectrodes=[505:508 497:499 4 8 489:492 7 484 11 19 478 14 38 470 10 42 13 379];
%StimElectrodes=[7 10 13 14 508];
StimElectrodes=[473 459 441 13];

StimTimes=[600:600:9600];

N=14
for i=1:N
    Movie=3+16*(i-1)
    [StimChannels,Amplitudes]=NS_StimulatedChannels('D:\Home\Pawel\analysis\2010-09-11-0\data003',13,Movie,[1:512],NS_GlobalConstants);
    A(i)=Amplitudes(1);
end
A=A/0.066*1.07;
    
AmplitudesPerFigure=3;
N0=ceil(N/AmplitudesPerFigure);
T=[1:Length]/20;
for SE=1:length(StimElectrodes);
    Info=NS512_MovieNumberAndTimingForPattern('G:\uncompressed\2010-09-14-0\movie002',StimElectrodes(SE));
    MovieStart=Info(1,1);
    
    XStart=Info(1,3);
    %T=([1:Length]-XStart)/20;
    if MovieStart==16
        N=N-1;
    end
    for M=1:N
        FigureNumber=ceil(M/AmplitudesPerFigure)
        SubplotNumber=M-(FigureNumber-1)*AmplitudesPerFigure
        
        Movie=MovieStart+16*(M-1)  
        f=fopen(['D:\Home\Pawel\analysis\2010-09-14-0\data002\mv_' num2str(Movie)]);
        a=fread(f,'int16');
        fclose(f);
        b=reshape(a,NumberOfElectrodes,NumberOfRepetitions,Length);
        s1=reshape(b(1,:,:),NumberOfRepetitions,Length);
        s2=reshape(b(2,:,:),NumberOfRepetitions,Length);    
        s3=reshape(b(3,:,:),NumberOfRepetitions,Length);    
        s4=reshape(b(4,:,:),NumberOfRepetitions,Length);    
        
        figure(FigureNumber);
    
        %hold on;    
        for i=StimTimes
            s1(:,i:i+10)=-300;
            s2(:,i:i+10)=-300;
            s3(:,i:i+10)=-300;
            s4(:,i:i+10)=-300;
            %h=plot([i i],[-250 -200],'r-');
        end
        subplot(3,AmplitudesPerFigure,SubplotNumber);
        plot(T-XStart/20,s1(1:50,:)');  
        grid on
        h=text(XStart+100,-250,num2str(Movie));
        set(h,'FontSize',20);
        axis([0 12 -400 -250]);
        xlabel('Time [ms]');
        ylabel('Output signal [mV]');
        text(9,-370,['I=' num2str(A(M),'%10.3f') '\muA']);
   
        sm=mean(s3,1);  
        for i=1:50
            s3(i,:)=s3(i,:)-sm;
        end
        subplot(3,AmplitudesPerFigure,SubplotNumber+AmplitudesPerFigure);
        h=plot(T-XStart/20,s3');
        grid on
        %set(h,'Color','b')
        axis([0 12 -40 30]);
        xlabel('Time [ms]');
        ylabel('Output signal [mV]');
        for i=1:50
            s3(i,:)=s3(i,:)-sm;
        end
        
        subplot(3,AmplitudesPerFigure,SubplotNumber+2*AmplitudesPerFigure);
        h=plot(T-XStart/20,sm);
        grid on
        %set(h,'Color','b')
        axis([0 12 -360 -260]);    
        xlabel('Time [ms]');
        ylabel('Output signal [mV]');
    end
    for f=1:FigureNumber
        figure(f);
        FullName=['D:\Home\Pawel\analysis\2010-09-14-0\figures\stim' num2str(StimElectrodes(SE)) '_' num2str(f)];            
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]); 
        print(h, '-dtiff', '-r120', FullName);
    end
end

ChunkData=NS_MovieData(DataFilename,Movie,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);