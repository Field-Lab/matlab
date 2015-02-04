% clear all

params.saveFigures = 0;
params.shiftLimSpikeMin = [2 35];
params.shiftStep = 0.25;
params.residLim  = [6 40]; %stim ends at 250us
modelType = 'prevArtifact';
stimSystem = '512'; %'512' or '61'
%patternNos = [39 457 410 289]; %289 is control
patternNos = [74 70]; %289 is control

movieInt = 0;

%LOOP over all experiment folders

%for i=2:29; 
for i=1; 
	%% filling elecRespInfo with details for creation of elecResp

	dataPathFolder = ['data' sprintf('%3.3d',i)];
	disp(['Currently Analysing: ' dataPathFolder])
	elecRespInfo.experimentName = '2014-11-20-0';
	elecRespInfo.dataPath       = ['/Volumes/Analysis/2014-11-20-0/' dataPathFolder '/'];  %Location of raw data chunks
	elecRespInfo.analysisPath   = '/Volumes/Analysis/2014-11-20-0/data000/';  %Location of vision output files

	% Check inputs
	%{
	expName = input(['Experiment name is ' ...
			elecRespInfo.experimentName '. \n OK? Enter Y or correct name: '],'s');
	if ~strcmpi(expName,'Y')
			elecRespInfo.experimentName = expName;
			disp(['experiment name set to: ' elecRespInfo.experimentName '\n']);
	end

	otherPath = input(['\nLocation of the data separated by patterns is ' ...
			elecRespInfo.dataPath '. \n OK? Enter Y or desired path: '],'s');
	if ~strcmpi(otherPath,'Y') %compare strings ignore case
		 if exist(otherPath,'dir')
				 elecRespInfo.dataPath = otherPath;
				 disp(['data location set to: ' elecRespInfo.dataPath '\n']);
		 else
				 disp('invalid data location. End script');
				 return; 
		 end
	end
			
	otherPath = input(['\nLocation of the Vision output files is ' ...
			elecRespInfo.analysisPath '. \n OK? Enter Y or desired path: '],'s');
	if ~strcmpi(otherPath,'Y') 
		 if exist(otherPath,'dir')
				 elecRespInfo.analysisPath = otherPath;
				 disp(['analysis path set to: ' elecRespInfo.dataPath]);
		 else
				 disp('invalid analysis path. End script');
				 return; 
		 end
	end
	%}

	i = find(elecRespInfo.analysisPath == filesep,2,'last');

	elecRespInfo.analysisBaseName = elecRespInfo.analysisPath(i(1)+1:i(2)-1)
	elecRespInfo.savePath         = elecRespInfo.dataPath;
	elecRespInfo.movieInt         =   movieInt;

	%elecRespInfo.mainNeuron =       input('Enter main neuron id: '); %108;
	elecRespInfo.mainNeuron = 922; 
	elecRespInfo.activeNeurons{1} = [];
	elecRespInfo.activeNeurons{2} = [];
	elecRespInfo.activeNeurons{3} = [];
	elecRespInfo.activeNeurons{4} = [];

	elecRespInfo.pElec    =      []; %input('Enter primary stimulating electrode (if known): '); %   [];  %LG: primary STIMULATING  electrode. from data015 on 2011-06-24 elecResp.stimInfo.pElec = 0. 
	 
	%elecRespInfo.mainElec =      input(['Enter main recording electrode for neuron ' num2str(elecRespInfo.mainNeuron) ': ']); %2;  %LG: main RECORDING electrode
	elecRespInfo.mainElec =  62; 
	elecRespInfo.otherElecs{1} = [];   % Other recording electodes. 
	elecRespInfo.otherElecs{2} = [];
	elecRespInfo.otherElecs{3} = [];
	elecRespInfo.otherElecs{4} = [];
	elecRespInfo.otherElecs{5} = [];


	elecRespInfo.artifactPath =   []; %'/Volumes/Analysis/2014-04-15-4/data004/';
	elecRespInfo.artMovieFirst =  [];
	elecRespInfo.artMovieLast =   [];
	elecRespInfo.artMovieInt =    [];
	elecRespInfo.sampleRate =     20000;

	elecRespInfo.autoMovie = true;
	elecRespInfo.externalEi = false;

	elecRespInfo.stimSystem = stimSystem; 

	for i = 1:length(patternNos)

			elecRespInfo.patternNo =  patternNos(i);
			%elecRespInfo.movieFirst = firstMovies(i);
			%elecRespInfo.movieLast =  lastMovies(i);
			%elecRespInfo.movieFirst = 1;
			%elecRespInfo.movieLast =  1;
			
			elecRespInfo.shortName =  [elecRespInfo.experimentName '_', num2str(elecRespInfo.mainNeuron), '_p' num2str(patternNos(i))];

			elecResp = createElecRespStruct(elecRespInfo);

			%temp = load(['/Volumes/Palace/Analysis/Lauren/2008-08-26-0/data008/elecResp_n257_p' num2str(patternNos(i)) '.mat']);
			%elecResp = temp.elecResp;
			
			elecRespName = ['elecResp_n' num2str(elecRespInfo.mainNeuron) '_p' num2str(patternNos(i)) '.mat'];
			save([elecRespInfo.dataPath filesep elecRespName], 'elecResp')
			

			for j = 1:length(elecResp.stimInfo.movieNos)
					elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), params,...
							'modelType', modelType);
							
					save([elecResp.names.data_path filesep elecRespName], 'elecResp')
					disp(['done analyzing movie ' num2str(elecResp.stimInfo.movieNos(j)) ', pattern ' num2str(patternNos(i))])
			end
	end
end


