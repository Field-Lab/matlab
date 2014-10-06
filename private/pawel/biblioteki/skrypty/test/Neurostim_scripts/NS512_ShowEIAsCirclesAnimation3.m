RawDataPath='G:\analysis\2010-09-14-0\data000'; %define path to raw data file
paramsFilePath='C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_07_18\VisionOutput\data000\data000.params';
neuronFilePath='C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_07_18\VisionOutput\data000\data000.neurons';

NeuronEI=NS512_LoadNeuronEIFromRawData(RawDataPath,paramsFilePath,neuronFilePath,6841,5900,400,0);

for i=1:512
    NeuronEI(i,:)=NeuronEI(i,:)-NeuronEI(i,400);
end

FigureName='C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_07_18\other_plots\ble.gif';

Neurons=[7 11 22 24];
figure(1)
for Frame=1:400
    clf
    Frame
    for i=1%1:length(NeuronWithManyStimulations)
        Events=find(neurons==NeuronWithManyStimulations(i));
        S1=sigmas(Events);
        DirectivityIndex(i)=length(find(S1<7))/length(S1);    
    
        %clf;
        hold on;
        for j=1:length(Events)
            h15=plot([XStim(Events(j)) XSeed(Events(j))],[YStim(Events(j)) YSeed(Events(j))],'b-');
            text(XStim(Events(j)),YStim(Events(j)),num2str(j));
            if sigmas(Events(j))<40 %7
                set(h15,'Color','r');
            end        
        end
        axis([-1000 1000 -500 500]);
    end
            
    h=NS512_ShowEIAsCircles(2*NeuronEI(:,Frame),500,[1:512],[],[-1000 1000],[-500 500]);
    pause(0.01)
    
    h=gcf;
    frame = getframe(h);
    im = frame2im(frame);
    [imind,map] = rgb2ind(im,256);
    if Frame == 1
        imwrite(imind,map,FigureName,'gif', 'Loopcount',1);
    else
        imwrite(imind,map,FigureName,'gif','WriteMode','append','DelayTime',0.1);
    end            
end