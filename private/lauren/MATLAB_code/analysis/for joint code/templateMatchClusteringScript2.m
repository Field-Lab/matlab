clear all

%addpath(genpath('/snle/home/lhruby/MATLAB_code/analysis/from Pawel/matlab_cp_2008_11_19/biblioteki/'))

params.saveFigures = 0;
params.shiftLimSpikeMin = [6 35];
params.shiftStep = 0.25;
params.residLim = [9 40];
modelType = 'prevArtifact';

patternNos = 768:909;

movieInt = 0;

%% filling elecRespInfo with details for creation of elecResp

elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data003/';
%elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data002/';

elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data000-cf-lh/';
%elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data000/';
elecRespInfo.analysisBaseName = 'data000-cf-lh';
elecRespInfo.savePath =         elecRespInfo.dataPath;

elecRespInfo.movieInt =   movieInt;

%elecRespInfo.mainNeuron =       151;
%elecRespInfo.mainNeuron =       560;
elecRespInfo.mainNeuron =       632;
elecRespInfo.activeNeurons{1} = [];
elecRespInfo.activeNeurons{2} = [];
elecRespInfo.activeNeurons{3} = [];
elecRespInfo.activeNeurons{4} = [];

elecRespInfo.pElec =         43;
 
elecRespInfo.mainElec =      43;
elecRespInfo.otherElecs{1} = [];
elecRespInfo.otherElecs{2} = [];
elecRespInfo.otherElecs{3} = [];
elecRespInfo.otherElecs{4} = [];
elecRespInfo.otherElecs{5} = [];

elecRespInfo.experimentName = '2011-07-14-7';
elecRespInfo.artifactPath =   '';
elecRespInfo.artMovieFirst =  [];
elecRespInfo.artMovieLast =   [];
elecRespInfo.artMovieInt =    [];
elecRespInfo.sampleRate =     20000;

elecRespInfo.autoMovie = true;
elecRespInfo.externalEi = false;

for i = 1:length(patternNos)

    elecRespInfo.patternNo =  patternNos(i);
    %elecRespInfo.movieFirst = firstMovies(i);
    %elecRespInfo.movieLast =  lastMovies(i);
    elecRespInfo.movieFirst = 1;
    elecRespInfo.movieLast =  1;
    
    elecRespInfo.shortName =  ['2011-07-14-7_data003_n', num2str(elecRespInfo.mainNeuron), '_p' num2str(patternNos(i))];

    elecResp = createElecRespStruct(elecRespInfo);

    elecRespName = ['elecResp_n' num2str(elecRespInfo.mainNeuron) '_p' num2str(patternNos(i)) '.mat'];
    save([elecRespInfo.dataPath filesep elecRespName], 'elecResp')
    

    for j = 1:length(elecResp.stimInfo.movieNos)
        elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), params,...
            'modelType', modelType);
            
        save([elecResp.names.data_path filesep elecRespName], 'elecResp')
        disp(['done analyzing movie ' num2str(elecResp.stimInfo.movieNos(j)) ', pattern ' num2str(patternNos(i))])
    end
end


