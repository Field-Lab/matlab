Movie=41;
NumberOfElectrodes=14;
NumberOfRepetitions=50;
Length=10000;
Fi=21
Ch1=20;
Ch2=446;
ChannelsPlot1=electrodeMap.getAdjacentsTo(Ch1,1)';
ChannelsPlot2=electrodeMap.getAdjacentsTo(Ch2,1)';
Channels=[ChannelsPlot1 ChannelsPlot2];
%21,1
for i=122%:8:135
        f=fopen(['D:\Home\Pawel\analysis\2010-09-14-0\data008new\mv_' num2str(i)]);
        a=fread(f,'int16');
        fclose(f);
        b=reshape(a,NumberOfElectrodes,NumberOfRepetitions,Length);
        figure(Fi)
        
        for j=1:7
            subplot(3,3,j);
            s1=reshape(b(j,:,:),NumberOfRepetitions,Length);        
            h=plot(s1');
            set(h,'Color','b');
            axis([1200 1700 -420 -260])
            grid on
        end
        
        %subplot(2,1,2);
        %s1=reshape(b(2,:,:),NumberOfRepetitions,Length);        
        %h=plot(s1');
        %set(h,'Color','b');
        % axis([0 5500 -340 -260])
        %grid on
        %pause(3);
end
break
FullName=['D:\Home\Pawel\analysis\2010-09-14-0\figures\data008'];            
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]); 
        print(h, '-dtiff', '-r120', FullName);