clear all

params.saveFigures = 0;
params.shiftLimSpikeMin = [3 35];
params.shiftStep = 0.25;
params.residLim = [7 40];
%modelType = 'ttx';
modelType = 'prevArtifact';

%nOccurrences = [76 88 100 104 88 90];
%occurancesAll = {[28:54 104:125 148]; [10:18 34:39 46]; [29:56 109:132 157]; [31:60 112:132 154];...
%    [40:78 155:191 229]; [31:60 120:148 178]; [27:52 102:124 148]; [27:52 104:128 154]};

occurancesAll = {[1:20 41:58]; [1:23 47:67]; [1:26 53:76]; [1:25 51:77]; [1:22 45:66]; [1:26 53:71]};
%occurancesAll = {[]; []; []; []; []; [145:148 178]; [27:36 51:52 102:123]; [44:52 104:114 121:124 128 154]};

for patternNo = [2:3]
    %occurances = 1:nOccurrences(patternNo);
    occurances = occurancesAll{patternNo};
    
    movieInt = 0;
    
    %% filling elecRespInfo with details for creation of elecResp
    
    % Points to preprocessed electrical stimulation data.
    elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2013-05-28-3/data012/'; 
    %elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/data013/'; 
    %elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2013-05-28-3/data012/';    
    %elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data015/';
    %elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/data013/';
    %elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data002/';
    
    % Points to white noise vision output.
    elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2013-05-28-3/data010/';
    %elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data006/';
    %elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/data005/';
    %elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data000/';
    elecRespInfo.analysisBaseName = 'data010';
    elecRespInfo.savePath =         elecRespInfo.dataPath;
    
    elecRespInfo.movieInt =   movieInt;
    
    %elecRespInfo.mainNeuron =       151;
    %elecRespInfo.mainNeuron =       560;
    %elecRespInfo.mainNeuron =       784;
    %elecRespInfo.mainNeuron =       498;
    %elecRespInfo.mainNeuron =       618;
    %elecRespInfo.mainNeuron =       918;
    %elecRespInfo.mainNeuron =       48;
    %elecRespInfo.mainNeuron =       187;
    %elecRespInfo.mainNeuron =       304;
    %elecRespInfo.mainNeuron =       407;
    elecRespInfo.mainNeuron =       166;
    
    elecRespInfo.activeNeurons{1} = [];
    elecRespInfo.activeNeurons{2} = [];
    elecRespInfo.activeNeurons{3} = [];
    elecRespInfo.activeNeurons{4} = [];
    
    elecRespInfo.pElec =         [];
    
    %elecRespInfo.mainElec =      50;
    %elecRespInfo.mainElec =      37;
    %elecRespInfo.mainElec =      45;
    %elecRespInfo.mainElec =      61;
    %elecRespInfo.mainElec =      4;
    %elecRespInfo.mainElec =      14;
    %elecRespInfo.mainElec =      24;
    %elecRespInfo.mainElec =      33;
    elecRespInfo.mainElec =      6;

    elecRespInfo.otherElecs{1} = [];
    elecRespInfo.otherElecs{2} = [];
    elecRespInfo.otherElecs{3} = [];
    elecRespInfo.otherElecs{4} = [];
    elecRespInfo.otherElecs{5} = [];
    
    %elecRespInfo.experimentName = '2011-07-14-0';
    elecRespInfo.experimentName = '2013-05-28-3';
    %elecRespInfo.artifactPath =   '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/data015/';
    elecRespInfo.artifactPath =   '/snle/lab/Experiments/Array/Analysis/2013-05-28-3/data015/';
    elecRespInfo.artMovieFirst =  [];
    elecRespInfo.artMovieLast =   [];
    elecRespInfo.artMovieInt =    [];
    elecRespInfo.sampleRate =     20000;
    
    elecRespInfo.autoMovie = true;
    elecRespInfo.externalEi = false;
    
    for i = 1:length(occurances)
        
        elecRespInfo.patternNo =  [num2str(patternNo) '_o' num2str(occurances(i))];
        elecRespInfo.movieFirst = 1;
        elecRespInfo.movieLast =  1;
        
        elecRespInfo.shortName =  ['2013-05-28-3_data012_n', num2str(elecRespInfo.mainNeuron),...
            '_p' elecRespInfo.patternNo];
        
        elecResp = createElecRespStruct(elecRespInfo);
        
        elecRespName = ['elecResp_n' num2str(elecRespInfo.mainNeuron) '_p' elecRespInfo.patternNo '.mat'];
        save([elecRespInfo.dataPath filesep elecRespName], 'elecResp')
        
        
        for j = 1:length(elecResp.stimInfo.movieNos)
            if j == 1
               elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), params,...
                   'modelType', 'ttx', 'spikeMinExclusionReg', 0);
            else
                elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), params,...
                    'modelType', modelType, 'spikeMinExclusionReg', 0);
                
                if elecResp.analysis.successRates(j) == 1
                    disp('success rate = 1; reanalyzing with ttx artifact')
                    elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), params,...
                        'modelType', 'ttx', 'spikeMinExclusionReg', 0);

                    disp('success rate = 1; reanalyzing')
%                     for kk = 1:10
%                         elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), params,...
%                             'modelType', 'currentArtifact', 'spikeMinExclusionReg', 0);
%                     end
                end
                
            end
            
            
            save([elecResp.names.data_path filesep elecRespName], 'elecResp')
            disp(['done analyzing movie ' num2str(elecResp.stimInfo.movieNos(j)) ', pattern ' elecRespInfo.patternNo])
        end
    end
end

%% 
