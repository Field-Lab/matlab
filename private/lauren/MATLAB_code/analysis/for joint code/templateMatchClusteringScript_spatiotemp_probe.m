clear all

params.saveFigures = 0;
params.shiftLimSpikeMin = [7 35];
params.shiftStep = 0.25;
params.residLim = [11 40];
modelType = 'prevArtifact';


primElec = 60;
prePulseElecs = [60];
delays = [10 20 40 80 160];
preAmps = {'0_703', '0_803', '2_008'}; %write as strings with _ replacing . (should match format of preprocessed data files)

movieInt = 0;

analyzePrepulse = false;

%% filling elecRespInfo with details for creation of elecResp

elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data004/';
%elecRespInfo.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data002/';

elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data000/';
%elecRespInfo.analysisPath =     '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data000/';
elecRespInfo.analysisBaseName = 'data000';
elecRespInfo.savePath =         elecRespInfo.dataPath;

elecRespInfo.movieInt =   movieInt;

%elecRespInfo.mainNeuron =       151;
%elecRespInfo.mainNeuron =       560;
elecRespInfo.mainNeuron =       948;
elecRespInfo.activeNeurons{1} = [];
elecRespInfo.activeNeurons{2} = [];
elecRespInfo.activeNeurons{3} = [];
elecRespInfo.activeNeurons{4} = [];

elecRespInfo.pElec =         60;
 
elecRespInfo.mainElec =      60;
elecRespInfo.otherElecs{1} = [];
elecRespInfo.otherElecs{2} = [];
elecRespInfo.otherElecs{3} = [];
elecRespInfo.otherElecs{4} = [];
elecRespInfo.otherElecs{5} = [];

elecRespInfo.experimentName = '2011-08-04-5';
elecRespInfo.artifactPath =   '';
elecRespInfo.artMovieFirst =  [];
elecRespInfo.artMovieLast =   [];
elecRespInfo.artMovieInt =    [];
elecRespInfo.sampleRate =     20000;

elecRespInfo.autoMovie = true;
elecRespInfo.externalEi = false;

elecRespInfo.movieFirst = 1;
elecRespInfo.movieLast =  1;

elecRespInfo.isSpatioTempProbe = true;

for zz = 1:2
    if zz == 1
        analyzePrepulse = false;
    else
        analyzePrepulse = true;
        delays(delays==0) = [];
    end
    
    for kk = 1:length(prePulseElecs)
        for ii = 1:length(delays)
            for jj = 1:length(preAmps)
                
                %shift analysis sample limits to where analyzed pulse is within
                %the trace
                paramsCopy = params;
                if ~analyzePrepulse
                    paramsCopy.shiftLimSpikeMin = params.shiftLimSpikeMin + delays(ii);
                    paramsCopy.residLim = params.residLim + delays(ii);
                elseif delays(ii) == 20
                    paramsCopy.shiftLimSpikeMin(2) = 20;
                    paramsCopy.residLim(2) = 22;
                end
                
%                 if delays(ii) == 10 && analyzePrepulse
%                     elecRespInfo.mainElec = 22;
%                 else
%                     elecRespInfo.mainElec = 17;
%                 end
                
                elecRespInfo.patternNo =  [num2str(primElec) '_pre' num2str(prePulseElecs(kk)) '_d' num2str(delays(ii)) '_a' preAmps{jj}];
                
                elecRespInfo.shortName =  ['2011-08-04-5_data004_n', num2str(elecRespInfo.mainNeuron),...
                    '_p' elecRespInfo.patternNo];
                
                elecResp = createElecRespStruct(elecRespInfo);
                
                if ~analyzePrepulse
                    elecRespName = ['elecResp_n' num2str(elecRespInfo.mainNeuron) '_p' elecRespInfo.patternNo '.mat'];
                else
                    elecRespName = ['elecResp_n' num2str(elecRespInfo.mainNeuron) '_p' elecRespInfo.patternNo '_preResponse.mat'];
                end
                
                save([elecRespInfo.dataPath filesep elecRespName], 'elecResp')
                
                
                for j = 1:length(elecResp.stimInfo.movieNos)
                    elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), paramsCopy,...
                        'modelType', modelType);
                    
                    save([elecResp.names.data_path filesep elecRespName], 'elecResp')
                    disp(['done analyzing movie ' num2str(elecResp.stimInfo.movieNos(j)) ', pattern ' elecRespInfo.patternNo])
                end
            end
        end
    end
end

