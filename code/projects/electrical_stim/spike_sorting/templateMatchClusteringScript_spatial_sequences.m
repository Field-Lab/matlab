clear all

params.saveFigures = 0;
params.shiftLimSpikeMin = [7 45];
params.shiftStep = 0.25;
params.residLim = [11 50];
modelType = 'prevArtifact';


patternNo = 14;
occurances = 1:4;
delays = [20 40 100];

movieInt = 0;

%% filling elecRespInfo with details for creation of elecResp
elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data010/';
%elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/data013/';
%elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data002/';

elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data009/';
%elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/data005/';
%elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data000/';
elecRespInfo.analysisBaseName = 'data009';
elecRespInfo.savePath =         elecRespInfo.dataPath;

elecRespInfo.movieInt =   movieInt;

%elecRespInfo.mainNeuron =       151;
%elecRespInfo.mainNeuron =       560;
elecRespInfo.mainNeuron =       228;
elecRespInfo.activeNeurons{1} = [];
elecRespInfo.activeNeurons{2} = [];
elecRespInfo.activeNeurons{3} = [];
elecRespInfo.activeNeurons{4} = [];

elecRespInfo.pElec =         14;
 
elecRespInfo.mainElec =      13;
elecRespInfo.otherElecs{1} = [];
elecRespInfo.otherElecs{2} = [];
elecRespInfo.otherElecs{3} = [];
elecRespInfo.otherElecs{4} = [];
elecRespInfo.otherElecs{5} = [];

elecRespInfo.experimentName = '2011-06-24-5';
elecRespInfo.artifactPath =   '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data010/';
elecRespInfo.artMovieFirst =  [];
elecRespInfo.artMovieLast =   [];
elecRespInfo.artMovieInt =    [];
elecRespInfo.sampleRate =     20000;

elecRespInfo.autoMovie = true;
elecRespInfo.externalEi = false;

elecRespInfo.movieFirst = 1;
elecRespInfo.movieLast =  1;


for i = 1:length(occurances)
    for kk = 1:length(delays)

    elecRespInfo.patternNo =  [num2str(patternNo) '_d' num2str(delays(kk)) '_o' num2str(occurances(i))];
    
    elecRespInfo.shortName =  ['2011-06-24-5_data010_n', num2str(elecRespInfo.mainNeuron),...
        '_p' elecRespInfo.patternNo];

    elecResp = createElecRespStruct(elecRespInfo);
    
    elecRespName = ['elecResp_n' num2str(elecRespInfo.mainNeuron) '_p' elecRespInfo.patternNo '.mat'];
    save([elecRespInfo.dataPath filesep elecRespName], 'elecResp')
    

    for j = 1:length(elecResp.stimInfo.movieNos)
        elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), params,...
            'modelType', modelType);
        
%         if elecResp.analysis.successRates(j) == 1
%             elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), params,...
%                 'modelType', 'ttx', 'noRefit', true);
%         else
%             elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), params,...
%                 'modelType', 'currentArtifact');
%         end
            
        save([elecResp.names.data_path filesep elecRespName], 'elecResp')
        disp(['done analyzing movie ' num2str(elecResp.stimInfo.movieNos(j)) ', pattern ' elecRespInfo.patternNo])
    end
    end
end


