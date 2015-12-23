clear all


%params.shiftLimSpikeMin = [6 35];
%params.residLim = [8 40];

elecRespInfo.patternNo = 6;
pw = 50;

%excludeNeurons = [496 512 541];
%moviesToAnalyze = 24:26;



%excludeNeurons = [91 286 302 391 424 556 586 631 767 856 902 48 62 111 153 ...
%    169 187 259 305 338 381 395 444 457 470 578 624 647 680 754 755 799 814 877 879 951 2 152 257 407 436 453 661 768 871 18 138 288 392 573 617];

%excludeNeurons = [91 286 302 391 424 556 586 631 767 856 902 48 62 111 153 ...
%    169 187 259 305 338 381 395 444 457 470 578 624 647 680 754 755 799 814 877 879 951];

excludeNeurons = [91 286 302 391 424 556 586 631 767 856 902 48 62 111 153 169 187 259 305 338 381 395 444 457 470 578 624 647 680 754 951];

%excludeNeurons = [91 286 302 391 424 556 586 631 767 856 902 48 62 111 153 169 187 259 305 ...
%    338 381 395 444 457 470 578 624 647 680 754 755 799 814 877 879 951 2 152 257 407 436 453 661 768 871 18 138 288 392 573 617 751 753];
moviesToAnalyze = [31 33 35];

                                    
%% stuff that doesn't usually change


if 1 % for 2011-01-11-0/data029 (on laptop)
    elecRespInfo.dataPath = '/Analysis/lauren/2011-01-11-0/data029/';
    elecRespInfo.analysisPath = '/Analysis/lauren/2011-01-11-0/data030-lh/';
    elecRespInfo.analysisBaseName = 'data030-lh';
    elecRespInfo.experimentName = '2011-01-11-0';
    shortName_prefix =  '2011-01-11-0_data029_n';
    
    for ii = 1 %for code folding purposes only
        neuronRecElecs = [91    7;
            286     20;
            302     21;
            391     27;
            424     29;
            556     38;
            586     40;
            631     43;
            767     52;
            856     58;
            902     61;
            48      58;
            62      2;
            111     8;
            153     11;
            169     12;
            187     13;
            259     18;
            305     21;
            338     23;
            381     26;
            395     27;
            444     27;
            457     31;
            470     32;
            578     39;
            624     42;
            647     47;
            680     46;
            754     51;
            755     51;
            799     54;
            814     55;
            877     59;
            879     59;
            951     64;
            2       1;
            152     11;
            257     18
            407     28;
            436     30;
            453     31;
            661     45;
            768     52;
            871     59;
            18      2;
            138     10;
            288     20;
            392     27;
            573     39;
            617     42;
            751     51;
            753     48;
            888     60];
    end
    
end


if 0 % for 2008-08-27-2/data002 (on laptop)
    elecRespInfo.dataPath = '/Analysis/lauren/2008-08-27-2/data002/'; %#ok<UNRCH>
    elecRespInfo.analysisPath = '/Analysis/lauren/2008-08-27-2/data001-lh/';
    elecRespInfo.analysisBaseName = 'data001-lh';
    elecRespInfo.experimentName = '2008-08-27-2';
    shortName_prefix =  '2008-08-27-2_data002_n';
    
    for ii = 1 %for code folding purposes only
        neuronRecElecs = [18    2;
            91    7;
            92    7;
            138   10;
            166   12;
            167   12;
            182   13;
            241   17;
            256   18;
            274   19;
            332   23;
            348   24;
            379   15;
            395   27;
            437   30;
            496   34;
            512   35;
            526   36;
            527   38;
            541   40;
            590   40;
            616   42;
            646   44;
            676   46;
            677   46;
            707   51;
            766   52;
            782   54;
            811   55;
            857   41;
            858   54;
            872   61;
            886   60;
            931   63;
            108   8;
            257   18;
            317   22;
            513   35;
            754   52;
            902   62];
    end
end


elecRespInfo.pElec = elecRespInfo.patternNo;

existingFiles = dir([elecRespInfo.dataPath '*.mat']);
              
%% independent of experiment
              
params.saveFigures = 0;
params.shiftStep = 0.25;
modelType = 'prevArtifact';

if pw == 50
    elecRespInfo.autoMovie = false;
    elecRespInfo.movieInt =   2;
    elecRespInfo.movieFirst = 1;
    elecRespInfo.movieLast =  65;
    
    params.shiftLimSpikeMin = [6 35];
    params.residLim = [8 40];
elseif pw == 100
    elecRespInfo.autoMovie = false;
    elecRespInfo.movieInt =   2;
    elecRespInfo.movieFirst = 2;
    elecRespInfo.movieLast =  64;
    
    params.shiftLimSpikeMin = [7 35];
    params.residLim = [11 40];
end


% elecRespInfo.movieInt =   0;
% elecRespInfo.movieFirst = 0;
% elecRespInfo.movieLast =  0;

elecRespInfo.activeNeurons{1} = [];
elecRespInfo.activeNeurons{2} = [];
elecRespInfo.activeNeurons{3} = [];
elecRespInfo.activeNeurons{4} = [];  

elecRespInfo.otherElecs{1} = [];
elecRespInfo.otherElecs{2} = [];
elecRespInfo.otherElecs{3} = [];
elecRespInfo.otherElecs{4} = [];
elecRespInfo.otherElecs{5} = [];

elecRespInfo.artifactPath =   '';
elecRespInfo.artMovieFirst =  [];
elecRespInfo.artMovieLast =   [];
elecRespInfo.artMovieInt =    [];

elecRespInfo.sampleRate = 20000;

elecRespInfo.savePath = elecRespInfo.dataPath;

%elecRespInfo.autoMovie = true;
elecRespInfo.externalEi = false;



%% filling elecRespInfo with details for creation of elecResp


for i = 1:size(neuronRecElecs,1);
    
    if ~any(excludeNeurons == neuronRecElecs(i,1))
        
        elecRespInfo.mainNeuron = neuronRecElecs(i,1); 
        elecRespInfo.mainElec = neuronRecElecs(i,2);

        elecRespInfo.shortName =  [shortName_prefix num2str(elecRespInfo.mainNeuron), '_p' num2str(elecRespInfo.patternNo) '_w' num2str(pw)];

        elecResp = createElecRespStruct(elecRespInfo);

        elecRespName = ['elecResp_n' num2str(elecRespInfo.mainNeuron) '_p' num2str(elecRespInfo.patternNo) '_w' num2str(pw) '.mat'];
        
        
        if exist([elecResp.names.data_path filesep elecRespName], 'file')
            keyboard
        else
            save([elecRespInfo.dataPath filesep elecRespName], 'elecResp')

            %check to make sure you're not accidentally overwriting another
            %elecResp file
        
            %for j = 1:length(elecResp.stimInfo.movieNos)
            for j = 1:length(moviesToAnalyze)
                elecResp = templateMatchClustering(elecResp, moviesToAnalyze(j), params,...
                    'modelType', modelType);

                save([elecResp.names.data_path filesep elecRespName], 'elecResp')
                disp(['done analyzing movie ' num2str(moviesToAnalyze(j)) ', pattern ' num2str(elecRespInfo.patternNo) '_w' num2str(pw)])
            end
        end
    end
end


