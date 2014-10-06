function cellInfo = cell_list_cluster_stim(badElecThresh)

x = 0;

%may actually have 200 trials for each pattern (only 100 analyzed?)
x = x+1;
cellInfo(x).id = 559;
cellInfo(x).type = 'onPar';
cellInfo(x).pElec = 44;
cellInfo(x).patternNos = [221:226 233:238 245:250 257:262 269:274 281:286 293:298 305:310 317:330]; %only patterns with negative primary
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data008';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data007-NW/axonDirData.mat';
cellInfo(x).withTriplets = false;
cellInfo(x).tripletsAnalyzed = false;
cellInfo(x).xPlotLim = [-1 1];
cellInfo(x).yPlotLim = [0 1];
cellInfo(x).pitch = 60;

%may actually have 200 trials for each pattern (only 100 analyzed?)
x = x+1;
cellInfo(x).id = 31;
cellInfo(x).type = 'onMidg';
cellInfo(x).pElec = 3;
cellInfo(x).patternNos = [1:6 13:18 25:30 37:42 49:54 61:66 73:78 85:90 97:110]; %only patterns with negative primary
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data008';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data007-NW/axonDirData.mat';
cellInfo(x).withTriplets = false;
cellInfo(x).tripletsAnalyzed = false;
cellInfo(x).xPlotLim = [-0.8 0.8];
cellInfo(x).yPlotLim = [0 0.8];
cellInfo(x).pitch = 60;


% actual number of trials = 200 for each pattern (for most patterns, only
% the first 100 are analyzed)
x = x+1;
cellInfo(x).id = 800;
cellInfo(x).type = 'offPar';
cellInfo(x).pElec = 49;
cellInfo(x).patternNos = [221:226 233:238 245:250 257:262 269:274 281:286 293:298 305:310 317:330]; %only patterns with negative primary
cellInfo(x).excludePatterns = []; %exclude from line-fitting to determine lambda values
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data005';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data004-NW/axonDirData.mat';
cellInfo(x).withTriplets = false;
cellInfo(x).tripletsAnalyzed = false;
cellInfo(x).xPlotLim = [-1 1];
cellInfo(x).yPlotLim = [0 1];
cellInfo(x).pitch = 60;

%may actually have 200 trials for each pattern (only 100 analyzed?)
x = x+1;
cellInfo(x).id = 241;
cellInfo(x).type = 'onPar';
cellInfo(x).pElec = 21;
cellInfo(x).patternNos = [1:6 13:18 25:30 37:42 49:54 61:66 73:78 85:90 97:110];
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data002';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data000/axonDirData.mat';
cellInfo(x).withTriplets = false;
cellInfo(x).tripletsAnalyzed = false;
cellInfo(x).xPlotLim = [-4 4];
cellInfo(x).yPlotLim = [1 5];
cellInfo(x).pitch = 30;


x = x+1;
cellInfo(x).id = 34;
cellInfo(x).type = 'onPar';
cellInfo(x).pElec = 58;
cellInfo(x).patternNos = 559:837;
cellInfo(x).excludePatterns = [];
%cellInfo(x).excludePatterns = [583 584 587 595 596 599 600];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data032';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data030/axonDirData.mat';
cellInfo(x).withTriplets = true;
cellInfo(x).tripletsAnalyzed = true;
cellInfo(x).xPlotLim = [-1 1];
cellInfo(x).yPlotLim = [0 1];
cellInfo(x).pitch = 60;

x = x+1;
cellInfo(x).id = 1;
cellInfo(x).type = 'offPar';
cellInfo(x).pElec = 58;
cellInfo(x).patternNos = 559:837;
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data032';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data030/axonDirData.mat';
cellInfo(x).withTriplets = true;
cellInfo(x).tripletsAnalyzed = true;
cellInfo(x).xPlotLim = [-1.5 1.5];
cellInfo(x).yPlotLim = [0 1.5];
cellInfo(x).pitch = 60;

%unfinished analysis
% x = x+1;
% cellInfo(x).id = 166;
% cellInfo(x).type = '';
% cellInfo(x).pElec = ;
% cellInfo(x).patternNos = [];
% cellInfo(x).excludePatterns = [];
% cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/';
% cellInfo(x).withTriplets = ;
% cellInfo(x).tripletsAnalyzed = ;
% cellInfo(x).xPlotLim = [];
% cellInfo(x).yPlotLim = [];
% cellInfo(x).pitch = ;
%unfinished analysis

%
% evidence for poor recording electrode injecting less than expected
% current
x = x+1;
cellInfo(x).id = 725;
cellInfo(x).type = 'onPar';
cellInfo(x).pElec = 46;
cellInfo(x).patternNos = 304:606;
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data003_004/';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data001/axonDirData.mat';
cellInfo(x).withTriplets = true;
cellInfo(x).tripletsAnalyzed = true;
cellInfo(x).xPlotLim = [-5 5];
cellInfo(x).yPlotLim = [0 5];
cellInfo(x).pitch = 30; 

x = x+1;
cellInfo(x).id = 212;
cellInfo(x).type = 'onPar';
cellInfo(x).pElec = 14;
cellInfo(x).patternNos = 1:303;
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data003_004/';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data001/axonDirData.mat';
cellInfo(x).withTriplets = true;
cellInfo(x).tripletsAnalyzed = true;
cellInfo(x).xPlotLim = [-2.5 2.5];
cellInfo(x).yPlotLim = [0 2.5];
cellInfo(x).pitch = 30; 


x = x+1;
cellInfo(x).id = 482;
cellInfo(x).type = 'offPar';
cellInfo(x).pElec = 47;
cellInfo(x).patternNos = 1:87;
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-07-05-5/data004';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-07-05-5/data000/axonDirData.mat';
cellInfo(x).withTriplets = false;
cellInfo(x).tripletsAnalyzed = false;
cellInfo(x).xPlotLim = [-2.5 2.5];
cellInfo(x).yPlotLim = [0 2.5];
cellInfo(x).pitch = 30;

%shitty piece?  (really bad visual responses, sparse activity)
x = x+1;
cellInfo(x).id = 528;
cellInfo(x).type = 'onMidg';
cellInfo(x).pElec = 41;
cellInfo(x).patternNos = 304:606;
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data003';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data000/axonDirData.mat';
cellInfo(x).withTriplets = true;
cellInfo(x).tripletsAnalyzed = true;
cellInfo(x).xPlotLim = [-2.5 2.5];
cellInfo(x).yPlotLim = [0 2.5];
cellInfo(x).pitch = 60;

%shitty piece?  (really bad visual responses, sparse activity)
x = x+1;
cellInfo(x).id = 632;
cellInfo(x).type = 'onMidg';
cellInfo(x).pElec = 43;
cellInfo(x).patternNos = [607:678 896:909]; %only includes pairs
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data003';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data000/axonDirData.mat';
cellInfo(x).withTriplets = true;
cellInfo(x).tripletsAnalyzed = false;
cellInfo(x).xPlotLim = [-2.5 2.5];
cellInfo(x).yPlotLim = [0 2.5];
cellInfo(x).pitch = 60;

%shitty piece?  (really bad visual responses, sparse activity)
x = x+1;
cellInfo(x).id = 711;
cellInfo(x).type = 'offPar';
cellInfo(x).pElec = 43;
cellInfo(x).patternNos = [607:678 896:909]; %only includes pairs
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data003';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-07-14-7/data000/axonDirData.mat';
cellInfo(x).withTriplets = true;
cellInfo(x).tripletsAnalyzed = false;
cellInfo(x).xPlotLim = [-7 4];
cellInfo(x).yPlotLim = [0 5.5];
cellInfo(x).pitch = 60;

x = x+1;
cellInfo(x).id = 79;
cellInfo(x).type = 'offMidg';
cellInfo(x).pElec = 12;
cellInfo(x).patternNos = 1:87;
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-08-01-0/data005';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-08-01-0/data000/axonDirData.mat';
cellInfo(x).withTriplets = false;
cellInfo(x).tripletsAnalyzed = false;
cellInfo(x).xPlotLim = [-12 8];
cellInfo(x).yPlotLim = [0 10];
cellInfo(x).pitch = 60;

x = x+1;
cellInfo(x).id = 571;
cellInfo(x).type = 'onPar';
cellInfo(x).pElec = 43;
cellInfo(x).patternNos = [1:72 290:303]; %only includes pairs
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data006_007';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data000/axonDirData.mat';
cellInfo(x).withTriplets = true;
cellInfo(x).tripletsAnalyzed = false;
cellInfo(x).xPlotLim = [-7 3];
cellInfo(x).yPlotLim = [0 5];
cellInfo(x).pitch = 30;

x = x+1;
cellInfo(x).id = 646;
cellInfo(x).type = 'onPar';
cellInfo(x).pElec = 39;
cellInfo(x).patternNos = 88:174;
cellInfo(x).excludePatterns = [];
cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-12-13-0/data009';
cellInfo(x).pathToAxonDirData = '/snle/lab/Experiments/Array/Analysis/2011-12-13-0/data005/axonDirData.mat';
cellInfo(x).withTriplets = false;
cellInfo(x).tripletsAnalyzed = false;
cellInfo(x).xPlotLim = [-4 4];
cellInfo(x).yPlotLim = [0 4];
cellInfo(x).pitch = 60;

%unfinished analysis
% x = x+1;
% cellInfo(x).id = 212;
% cellInfo(x).type = 'onPar';
% cellInfo(x).pElec = 14;
% cellInfo(x).patternNos = 1:303;
% cellInfo(x).excludePatterns = [];
% cellInfo(x).relAmps = [-1.5 -1 -0.5 -0.25 -0.1 0.1 0.25 0.5 1 1.5];
% cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data003_004';
% cellInfo(x).withTriplets = true;
% cellInfo(x).tripletsAnalyzed = false;
% cellInfo(x).xPlotLim = [];
% cellInfo(x).yPlotLim = [];
% cellInfo(x).pitch = 30;

%% determine which electrodes are potentially "bad"

dsInfo = badElecCheck('badElecThresh', badElecThresh, 'makePlots', true);
dspath = cell(1, length(dsInfo));
for ii = 1:length(dsInfo);
    z = strfind(dsInfo(ii).eiPath, '/Analysis/');
    dspath{ii} = dsInfo(ii).eiPath(1:z+22);
end


for ii = 1:length(cellInfo)
    dsMatch = cellfun(@(x) ~isempty(strfind(cellInfo(ii).pathToData, x)), dspath);
    if ~sum(dsMatch)==1
        error(['could not find an unambiguous dataset match for ' cellInfo(ii).pathTOData])
    end
    
    badElecs = dsInfo(dsMatch).badElecs;
    sElecs = getCluster(cellInfo(ii).pElec);
    sElecs = sElecs(2:end);
    
    cellInfo(ii).badSElecs = intersect(badElecs, sElecs);
end


%%


%template

% x = x+1;
% cellInfo(x).id = ;
% cellInfo(x).type = '';
% cellInfo(x).pElec = ;
% cellInfo(x).patternNos = [];
% cellInfo(x).excludePatterns = [];
% cellInfo(x).relAmps = [];
% cellInfo(x).pathToData = '/snle/lab/Experiments/Array/Analysis/';
% cellInfo(x).withTriplets = ;
% cellInfo(x).tripletsAnalyzed = ;
% cellInfo(x).xPlotLim = [];
% cellInfo(x).yPlotLim = [];
% cellInfo(x).pitch = ;
