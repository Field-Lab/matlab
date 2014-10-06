function [electrodes,Array,MovieChunks]=NS512_1el_scan_simple(TimeShiftInMs,DelayInMs,NumberOfSamples);
%"simple" means - we just define one 8-electrode cluster (one electrode for
%each chip), and do not do the "turn off one electrode of those 8" thing.
%There is also only one movie here, with 64 patterns.

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;

Array=zeros(512,64);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
x=zeros(1,512);
y=zeros(1,512);
for i=1:512
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);
end

dx=30;
dy=60; %this in fact should be different, but that's how it is in the electrode map

xmin=min(min(X));
ymin=min(min(Y));

c=find(X==xmin & Y==ymin); %this is the electrode in "lower up" corner ofthe array. We start the electrode jumpiong from here.

%First, find order of electrodes within the single 64-electrode cluster. We
%further divide it to four 16-electrode clusters (4x4). Now, two consecutive different rows
%are shifted by 30 microns (half the pitch).

% relative coordinates within 4x4 cluster:
x=[0 60 120 180 30 90 150 210 0 60 120 180 30 90 150 210]; %relative coordinates within 4x4 cluster
y=[0 0 0 0 60 60 60 60 120 120 120 120 180 180 180 180];

% relative shifts for each of the 4x4 cluster within the 8x8 cluster:
xshifts4x4=[0 240 0 240];
yshifts4x4=[0 0 240 240];

% shifts for eight 8x8 clusters:
xshifts8x8=[0 480 960 1440 0 480 960 1440];
yshifts8x8=[0 0 0 0 480 480 480 480];

%The indexes i and j are just for different electrodes within 4x4 cluster.
%Four such clusters form 8x8 cluster, and there are eight 8x8 clusters on
%the array :)
for row=1:4 % for each of row of the 4x4 cluster... 
    for column=1:4 % for each column of the 4x4 cluster...
        index=(row-1)*4+column % which electrode within the 4x4 cluster
        X0coord=xmin+x(index);
        Y0coord=ymin+y(index); % coordinates within 4x4 cluster
        for cluster4x4=1:4 % for each 4x4 cluster within the 8x8 cluster...                                                
            X1coord=X0coord+xshifts4x4(cluster4x4);
            Y1coord=Y0coord+yshifts4x4(cluster4x4); %coordinates within 8x8 cluster
            electrodes2=[];
            for cluster8x8=1:8
                xcoord=X1coord+xshifts8x8(cluster8x8);
                ycoord=Y1coord+yshifts8x8(cluster8x8);
                
                electrode=find(X==xcoord & Y==ycoord);
                electrodes2=[electrodes2 electrode];
            end
            electrodes2;
            %so now the "electrodes2" array includes 8 electrodes - one from each 8x8 cluster.
            %Excluding one of them, we get 8 possible patterns. Each will
            %be in separate movie.
            
            PatternIndex=index+(cluster4x4-1)*16;
            Array(electrodes2,PatternIndex)=1;
            
            %for movie=1:8
            %    PatternIndex=index+(cluster4x4-1)*16+(movie-1)*64;
            %    ElectrodesForMovie=electrodes2([1:movie-1 movie+1:8]);
            %    Array(ElectrodesForMovie,PatternIndex)=1;
            %end                                              
        end                                                 
    end
end

electrodes=[1:512];
%testing
for i=1:0
    a=find(Array(:,i)==1);
    if length(a)~=7
        i
        warning('Number of patterns with active electrode different than 7');
    else
        n='OK';
    end
end

MovieChunks=[1];
for i=1:1 %generate different movie chunks
    Patterns=[1:64]+(i-1)*64;
    Times=[TimeShift:Delay:TimeShift+Delay*63];
    Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
    MovieChunks=[MovieChunks Chunk];
end