ChipAddresses=[24:31];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

%cd F:\2010-07-29-0;
%FileName='002';
%WritePath='C:\home\pawel\2010\analysis\07_2010_Cultures\2010-07-29-0\files_new';

%WritePathFigs=WritePath;

FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-50 50],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
%[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants,0,[1:90],[1:5],[1:512]);

for i=3:8
    MoviesToRead=[i:8:i+27*8];
    for j=0:7
        PatternsToRead=[1:8]+j*8;
        
        cd F:\2010-07-29-0;
        FileName='002';
        WritePath='I:\analysis\2010-07-29-0\data';
        WritePathFigs=WritePath;
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v2(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants,0,MoviesToRead,PatternsToRead,[1:512]);
        
        cd G:\2010-07-29-0;
        FileName='018';
        WritePath='I:\analysis\2010-07-29-0\artifacts';
        WritePathFigs=WritePath;
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v2(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants,0,MoviesToRead,PatternsToRead,[1:512]);                
    end
end

break;

cd E:\pawel\data\in_vivo\2009-07-25\analysis\;
FileName='006';
WritePath='E:\pawel\analysis\in-vivo\2009-07-25\data006';
WritePathFigs=WritePath;

FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-50 50],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
el=[20 70 76 106 139 354 358 383 421 444 445 463];
%el=[16 25 26 36 65 81 100 105 120 144];
electrodes=[];
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
Radius=1;
for i=el
    Channels=electrodeMap.getAdjacentsTo(i,Radius)';
    electrodes=[electrodes Channels];
end
electrodes=unique(electrodes);

electrodes1=[3 17 25 70 74 75 79 81 83 84 86 90 92 95 96 98 99 114 151 157];
electrodes2=[158 160 199 200 206 212 215 217 221 224 227 233 236 243 246 248 252 254];
electrodes3=[263 265 267 269 270 273 277 287 289 290 293 301 305 310 314 315 319 320 327];
electrodes4=[330 337 344 348 350 358 360 362 366 376 378 381 382 384 385 386 387 393 402 420 421 430 459 466];
electrodes=[electrodes1 electrodes2 electrodes3 electrodes4];

electrodes=[370 373 376 385 393 396 401 403 405 410 413 414 416 419 420 421 429 431 446 452 460 461 467 471 475 490 493 496 500 501 506 512];

NS_GlobalConstants=NS_GenerateGlobalConstants(512);
%break
electrodes=424:468;
WritePath='D:\analysis\2009-11-20-0\data003c';
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants,0,[1:2],[1:5],electrodes);
break;
electrodes=469:512;
WritePath='D:\analysis\2009-11-20-0\data003d';
%[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants,0,[1:2],[1:5],electrodes);

PrepDataPath='D:\analysis\retina\2009-11-27-0\data001';
ClusterFileName='D:\analysis\retina\2009-11-27-0\data001\ClusterFile_001';
GoodChannels=[385:512];
Movies=[1:2:151];

Patterns=[1:64];
Responses=zeros(length(Patterns),512);
Noise=Responses;
tic
for PatternNumber=Patterns
    PatternNumber
    [a,b]=NS512_ApplyLinearArtifactModel2(PrepDataPath,PatternNumber,Movies,GoodChannels,0,0,ClusterFileName,NS_GlobalConstants);
    Responses(PatternNumber,:)=a;
    Noise(PatternNumber,:)=b;
end
toc
save Responses;
Save Noise;