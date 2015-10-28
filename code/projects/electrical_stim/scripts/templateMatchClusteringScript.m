% clear all

params.saveFigures = 0;
params.shiftLimSpikeMin = [2 40];
params.shiftStep = 0.25;
params.residLim  = [6 45]; %stim ends at 250us
modelType = 'prevArtifact';
stimSystem = '512'; %'512' or '61'
% patternNos = 1:303;

% patternNos = [246 237:239 228:232 220:224 212:216 204:208 196:200 188:192 180:184 ...
%    172:176 164:168 156:160 148:152 140:144 132:136]; 
% patternNos = [244:246 234:238 227:231 218:222 211:215 203:206 195:198 ...
%     187:189 178:181 169:172 162:164 153:155 146:148 137 139 129:131]; 
patternNos = input('Enter pattern(s) to create elecResp file(s): ');

movieInt = 0;

%% filling elecRespInfo with details for creation of elecResp

% Experiment specific inputs
elecRespInfo.experimentName = '2015-09-23-3';
elecRespInfo.dataPath       = '/Volumes/Analysis/2015-09-23-3/data001/';  %Location of raw data chunks
elecRespInfo.analysisPath   = '/Volumes/Analysis/2015-09-23-3/data005/';  %Location of vision output files

% elecRespInfo.experimentName = '2015-09-18-0';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2015-09-18-0/data001/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2015-09-18-0/data000/';  %Location of vision output files

% elecRespInfo.experimentName = '2014-11-24-2';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2014-11-24-2/data002/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2014-11-24-2/data008/';  %Location of vision output files



% elecRespInfo.experimentName = '2015-04-14-0';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2015-04-14-0/data002/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2015-04-14-0/data005/';  %Location of vision output files

% elecRespInfo.experimentName = '2014-11-24-2';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2014-11-24-2/data006/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2014-11-24-2/data008/';  %Location of vision output files

% elecRespInfo.experimentName = '2014-09-10-0';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2014-09-10-0/data006/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2014-09-10-0/data000/';  %Location of vision output files

% elecRespInfo.experimentName = '2012-09-24-3';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2012-09-24-3/data006/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2012-09-24-3/data007/';  %Location of vision output files

% elecRespInfo.experimentName = '2014-08-20-1';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2014-08-20-1/data003/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2014-08-20-1/data006/';  %Location of vision output files

% elecRespInfo.experimentName = '2012-09-24-0';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2012-09-24-0/data007/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2012-09-24-0/data006/';  %Location of vision output files

% elecRespInfo.experimentName = '2014-07-08-3';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2014-07-08-3/data005/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2014-07-08-3/data001/';  %Location of vision output files

% elecRespInfo.experimentName = '2012-09-24-3';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2012-09-24-3/data008/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2012-09-24-3/data007/';  %Location of vision output files

% elecRespInfo.experimentName = '2014-04-15-4';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2014-04-15-4/data003/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2014-04-15-4/data000/';  %Location of vision output files

% elecRespInfo.experimentName = '2014-04-10-0';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2014-04-10-0/data005/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2014-04-10-0/data003/';  %Location of vision output files

% elecRespInfo.experimentName = '2014-06-04-3';
% elecRespInfo.dataPath       = '/Volumes/Analysis/2014-06-04-3/data008/';  %Location of raw data chunks
% elecRespInfo.analysisPath   = '/Volumes/Analysis/2014-06-04-3/data006/';  %Location of vision output files


% Check inputs
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

i = find(elecRespInfo.analysisPath == filesep,2,'last');

elecRespInfo.analysisBaseName = elecRespInfo.analysisPath(i(1)+1:i(2)-1);
elecRespInfo.savePath         =  elecRespInfo.dataPath; 
elecRespInfo.movieInt         =   movieInt;

elecRespInfo.mainNeuron =       input('Enter main neuron id: '); %108;
elecRespInfo.activeNeurons{1} = [];
elecRespInfo.activeNeurons{2} = [];
elecRespInfo.activeNeurons{3} = [];
elecRespInfo.activeNeurons{4} = [];

elecRespInfo.pElec    =      []; %input('Enter primary stimulating electrode (if known): '); %   [];  %LG: primary STIMULATING  electrode. from data015 on 2011-06-24 elecResp.stimInfo.pElec = 0. 
 
elecRespInfo.mainElec =      input(['Enter main recording electrode for neuron ' num2str(elecRespInfo.mainNeuron) ': ']); %2;  %LG: main RECORDING electrode
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

disp(elecRespInfo); 

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
    save([elecRespInfo.savePath filesep elecRespName], 'elecResp')
    

    for j = 1:length(elecResp.stimInfo.movieNos)
        elecResp = templateMatchClustering(elecResp, elecResp.stimInfo.movieNos(j), params,...
            'modelType', modelType);
            
        save([elecResp.names.savePath filesep elecRespName], 'elecResp')
        disp(['done analyzing movie ' num2str(elecResp.stimInfo.movieNos(j)) ', pattern ' num2str(patternNos(i))])
    end
end


