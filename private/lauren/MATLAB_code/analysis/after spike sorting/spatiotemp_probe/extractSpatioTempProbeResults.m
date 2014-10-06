function prePulse = extractSpatioTempProbeResults(pathToElecResp, stimElec, neuronID, prePulse, offsets, varargin)


p = inputParser;

p.addRequired('pathToElecResp', @ischar)
p.addRequired('stimElec', @isnumeric)
p.addRequired('neuronID', @isnumeric)
p.addRequired('prePulse', @isstruct)
p.addRequired('offsets', @isnumeric)


p.addParamValue('minDataForFitting', 50, @isnumeric) %minimum number of total data points (across all amplitudes) required for erf to be fit
p.addParamValue('plotErfFits', false, @islogical)

%if no response probabilities dip below upper limit OR no response
%probabilities get above lower limit, no erf is fit
p.addParamValue('probLimsForFitting', [0.25 0.75], @isnumeric)

p.parse(pathToElecResp, stimElec, neuronID, prePulse, offsets, varargin{:})

minData = p.Results.minDataForFitting;
probLims = p.Results.probLimsForFitting;
plotErfFits = p.Results.plotErfFits;


%%

% simultaneous stimulation not included in plots
offsets(offsets == 0) = [];

nPreElecs = length(prePulse);
nOffsets = length(offsets);

%% cases with some successes to prepulse

%extract analysis information
for ii = 1:nPreElecs
    %ii = 1;
    
    disp([10, '******extracting data for electrode ', num2str(prePulse(ii).elec), '******', 10])
    
    prePulse(ii).preResp = cell(nOffsets, length(prePulse(ii).amps));
    prePulse(ii).stimResp = cell(nOffsets, length(prePulse(ii).amps));
    prePulse(ii).stimAmps = cell(nOffsets, length(prePulse(ii).amps));
    prePulse(ii).meanPreRespProb = zeros(1, length(prePulse(ii).amps));
    for kk = 1:length(prePulse(ii).amps)
        nPrePulses = 0;
        for jj = 1:nOffsets
            aString = num2str(prePulse(ii).amps(kk), '%0-4.3f');
            aString(strfind(aString,'.')) = '_'; %replace decimal with underscore to match the file name convention
            
            %extract information from prePulses and store in prePulse
            %struct
            if prePulse(ii).analyzedPreResponse(kk)==1
                load([pathToElecResp filesep 'elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                    num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString '_preResponse.mat']);
                
                nMovies = length(elecResp.stimInfo.movieNos);
                prePulse(ii).preResp{jj,kk} = cell(1,nMovies);
                
                for mm = 1:nMovies
                    prePulse(ii).preResp{jj,kk}{mm} = elecResp.analysis.latencies{mm}~=0;
                    
                    %all prePulse movies should be analyzed (required)
                    if isempty(elecResp.analysis.type{mm})
                        error([pathToElecResp filesep 'elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                            num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString '_preResponse.mat movie'...
                            num2str(elecResp.stimInfo.movieNos(mm)) ' hasn''t been analyzed - aborting'])
                    end
                    %all prepulses should be finalized!
                    if ~elecResp.analysis.finalized(mm)
                        warning([pathToElecResp filesep 'elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                            num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString '_preResponse.mat movie'...
                            num2str(elecResp.stimInfo.movieNos(mm)) ' hasn''t been finalized']) %#ok<WNTAG>
                    end
                    
                    %to determine mean response probability for each prepulse amplitude
                    prePulse(ii).meanPreRespProb(kk) = prePulse(ii).meanPreRespProb(kk) + sum(prePulse(ii).preResp{jj,kk}{mm});
                    nPrePulses = nPrePulses + length(prePulse(ii).preResp{jj,kk}{mm});
                end
                clear elecResp
            end
            
            %extract stim pulse analysis
            load([pathToElecResp filesep 'elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString '.mat']);
            nMovies = length(elecResp.stimInfo.movieNos);
            
            prePulse(ii).stimResp{jj,kk} = cell(1,nMovies);
            for mm = 1:nMovies
                if ~isempty(elecResp.analysis.type{mm})
                    if ~elecResp.analysis.finalized(mm)
                        warning([pathToElecResp filesep 'elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                            num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString '.mat movie '...
                            num2str(elecResp.stimInfo.movieNos(mm)) ' hasn''t been finalized']) %#ok<WNTAG>
                    end
                    prePulse(ii).stimResp{jj,kk}{mm} = elecResp.analysis.latencies{mm}~=0;
                end
            end
            prePulse(ii).stimAmps{jj,kk} = elecResp.stimInfo.stimAmps;
            
            clear elecResp
            
            %split stim data into groups: responds to prepulse vs. doesn't
            %respond to prepulse and fit curve
            nMovies = length(prePulse(ii).stimAmps{jj,kk});
            
            dataPostSucc = [];
            dataPostFail = [];
            
            % split data into groups and fit erfs
            mIndSucc = 0;
            mIndFail = 0;
            for mm = 1:nMovies
                if ~isempty(prePulse(ii).stimResp{jj,kk}{mm})
                    if prePulse(ii).analyzedPreResponse(kk)==1
                        preResp = prePulse(ii).preResp{jj,kk}{mm};
                    elseif prePulse(ii).analyzedPreResponse(kk)==2
                        preResp = true(size(prePulse(ii).stimResp{jj,kk}{mm}));
                    else
                        keyboard
                    end
                    if sum(preResp) ~= 0
                        mIndSucc = mIndSucc+1;
                        successesPostSucc = prePulse(ii).stimResp{jj,kk}{mm}(preResp);
                        elecRespPostSucc.analysis.latencies{mIndSucc} = successesPostSucc; %not actual latencies; only used in bootstrapping
                        dataPostSucc = [dataPostSucc zeros(3,1)]; %#ok<AGROW>
                        
                        dataPostSucc(1,end) = abs(prePulse(ii).stimAmps{jj,kk}(mm)); %amplitude
                        dataPostSucc(2,end) = sum(successesPostSucc)/length(successesPostSucc); %proportion of trials that are successes
                        dataPostSucc(3,end) = length(successesPostSucc); %number of trials in group
                    end
                    if sum(~preResp) ~= 0
                        mIndFail = mIndFail+1;
                        successesPostFail = prePulse(ii).stimResp{jj,kk}{mm}(~preResp);
                        elecRespPostFail.analysis.latencies{mIndFail} = successesPostFail; %not actual latencies; only used in bootstrapping
                        dataPostFail = [dataPostFail zeros(3,1)]; %#ok<AGROW>
                        
                        dataPostFail(1,end) = abs(prePulse(ii).stimAmps{jj,kk}(mm));
                        dataPostFail(2,end) = sum(successesPostFail)/length(successesPostFail);
                        dataPostFail(3,end) = length(successesPostFail);
                    end
                end
            end
            
            %post-success analysis
            if ~isempty(dataPostSucc) && any(dataPostSucc(2,:)>probLims(1)) && any(dataPostSucc(2,:)<probLims(2)) && sum(dataPostSucc(3,:))>=minData
                erfParams = erfFitter(dataPostSucc, 2, -1, 'makePlot', plotErfFits);
                prePulse(ii).erfParamsPostSucc{jj,kk} = erfParams;
                prePulse(ii).threshPostSucc(jj,kk) = -erfParams(2)/erfParams(1);
                clear erfParams
                
                elecRespPostSucc.stimInfo.movieNos = 1:size(dataPostSucc,2); %only used to determine number of movies
                for mm = 1:length(elecRespPostSucc.stimInfo.movieNos)
                    elecRespPostSucc.analysis.type{mm} = 'whatevs'; %only checked for emptiness
                end
                
                elecRespPostSucc.stimInfo.stimAmps = dataPostSucc(1,:);
                
                prePulse(ii).threshStdPostSucc(jj,kk) = bootstrapThresh(elecRespPostSucc, 100);
                
                if plotErfFits
                    axes(gca) %#ok<LAXES>
                    hold on
                    plot([prePulse(ii).threshPostSucc(jj,kk)+prePulse(ii).threshStdPostSucc(jj,kk)...
                        prePulse(ii).threshPostSucc(jj,kk)-prePulse(ii).threshStdPostSucc(jj,kk)], [0.5 0.5], 'k-')
                    hold off
                end
            else
                prePulse(ii).erfParamsPostSucc{jj,kk} = [NaN NaN];
                prePulse(ii).threshStdPostSucc(jj,kk) = NaN;
                if isempty(dataPostSucc) || sum(dataPostSucc(3,:)) < minData
                    prePulse(ii).threshPostSucc(jj,kk) = NaN;
                    disp(['***elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                        num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString ' post-success '...
                        'had too few data points (< ' num2str(minData) ') to fit erf'])
                elseif ~any(dataPostSucc(2,:) > probLims(1))
                    prePulse(ii).threshPostSucc(jj,kk) = inf;
                    if offsets(jj) ~= 10 %this is expected to occur for 10 sample delay so don't display warning
                        disp(['***elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                            num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString ' post-success '...
                            'had no response probabilities above ' num2str(probLims(1))])
                    end
                elseif ~any(dataPostSucc(2,:) < probLims(2))
                    prePulse(ii).threshPostSucc(jj,kk) = -inf;
                    disp(['***elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                        num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString ' post-success '...
                        'had no response probabilities below ' num2str(probLims(2))])
                else error('There''s a bug in here somewhere...')
                end
            end
            
            %warns if success rate of post-success, 0.5 ms delay stim pulse
            %reaches above the lower response probability limit
            if offsets(jj) == 10 && ~isempty(dataPostSucc) && any(dataPostSucc(2,:) > probLims(1))
                disp(['***elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                    num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString...
                    ' has at least one post-success response probability above ' num2str(probLims(1))])
            end
            
            % post-failure analysis
            if prePulse(ii).analyzedPreResponse(kk)==1
                if any(dataPostFail(2,:)>probLims(1)) && any(dataPostFail(2,:)<probLims(2)) && sum(dataPostFail(3,:)) >= minData
                    erfParams = erfFitter(dataPostFail, 2, -1, 'makePlot', plotErfFits);
                    prePulse(ii).erfParamsPostFail{jj,kk} = erfParams;
                    prePulse(ii).threshPostFail(jj,kk) = -erfParams(2)/erfParams(1);
                    clear erfParams
                    
                    elecRespPostFail.stimInfo.movieNos = 1:size(dataPostFail,2); %only used to determine number of movies
                    for mm = 1:length(elecRespPostFail.stimInfo.movieNos)
                        elecRespPostFail.analysis.type{mm} = 'whatevs'; %only checked for emptiness
                    end
                    
                    elecRespPostFail.stimInfo.stimAmps = dataPostFail(1,:);
                    
                    prePulse(ii).threshStdPostFail(jj,kk) = bootstrapThresh(elecRespPostFail, 100);
                    
                    if plotErfFits
                        axes(gca) %#ok<LAXES>
                        hold on
                        plot([prePulse(ii).threshPostFail(jj,kk)+prePulse(ii).threshStdPostFail(jj,kk)...
                            prePulse(ii).threshPostFail(jj,kk)-prePulse(ii).threshStdPostFail(jj,kk)], [0.5 0.5], 'k-')
                        hold off
                    end
                else
                    prePulse(ii).erfParamsPostFail{jj,kk} = [NaN NaN];
                    prePulse(ii).threshStdPostFail(jj,kk) = NaN;
                    if sum(dataPostFail(3,:)) < minData
                        prePulse(ii).threshPostFail(jj,kk) = NaN;
                        disp(['***elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                            num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString ' post-failure '...
                            'had too few data points (< ' num2str(minData) ') to fit erf'])
                    elseif ~any(dataPostFail(2,:) > probLims(1))
                        prePulse(ii).threshPostFail(jj,kk) = inf;
                        disp(['***elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                            num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString ' post-failure '...
                            'had no response probabilities above ' num2str(probLims(1))])
                    elseif ~any(dataPostFail(2,:) < probLims(2))
                        prePulse(ii).threshPostFail(jj,kk) = -inf;
                        disp(['***elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                            num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString ' post-failure '...
                            'had no response probabilities below ' num2str(probLims(2))])
                    else error('There''s a bug in here somewhere...')
                    end
                end
            else
                disp(['no data-splitting performed on elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '_pre'...
                    num2str(prePulse(ii).elec) '_d' num2str(offsets(jj)) '_a' aString ': responses to prePulses are essentially all successes'])
                prePulse(ii).erfParamsPostFail{jj,kk} = [];
                prePulse(ii).threshPostFail(jj,kk) = 0;
                prePulse(ii).threshStdPostFail(jj,kk) = 0;
            end
        end
        if prePulse(ii).analyzedPreResponse(kk)==1
            prePulse(ii).meanPreRespProb(kk) = prePulse(ii).meanPreRespProb(kk)/nPrePulses;
        elseif prePulse(ii).analyzedPreResponse(kk)==2
            prePulse(ii).meanPreRespProb(kk) = 1;
        end
    end
end


%%
% for ii = 1:nPreElecs
%     prePulse(ii).meanPreRespProb = zeros(1, length(prePulse(ii).amps));
%     for kk = 1:length(prePulse(ii).amps)
%         if prePulse(ii).analyzedPreResponse(kk)==1
%             nPrePulses = 0;
%             for jj = 1:nOffsets
%                 %to determine mean response probability for each prepulse amplitude
%                 for mm = 1:length(prePulse(ii).preResp{jj,kk})
%                     prePulse(ii).meanPreRespProb(kk) = prePulse(ii).meanPreRespProb(kk) + sum(prePulse(ii).preResp{jj,kk}{mm});
%                     nPrePulses = nPrePulses + length(prePulse(ii).preResp{jj,kk}{mm});
%                 end
%             end
%             prePulse(ii).meanPreRespProb(kk) = prePulse(ii).meanPreRespProb(kk)/nPrePulses;
%         else
%             prePulse(ii).meanPreRespProb(kk) = 1;
%         end
%     end
% end

    
% %% plotting preparations
% offsetsMs = offsets/20;
% markerSize = 4;
% 
% nPreElecs = length(prePulse);
% 
% offsets(offsets == 0) = [];
% 
% nPreElecs = length(prePulse);
% nOffsets = length(offsets);
% 
% plotColors(1,:) = [90 156 0]/255; %pale grass
% plotColors(2,:) = [255 124 59]/255; %salmon
% plotColors(3,:) = [101 52 255]/255; %purple
% plotColors(4,:) = [52 198 247]/255; %aqua
% plotColors(5,:) = [238 55 128]/255; %calm magenta
% 
% %load target cell ei from one of the elecResps
% temp = load([pathToElecResp filesep 'elecResp_n' num2str(neuronID) '_p' num2str(stimElec) '.mat']);
% threshIso = temp.elecResp.analysis.threshold;
% threshStdIso = temp.elecResp.analysis.threshStd;
% 
% for ii = 1:64
%     if ~(ii==9||ii==25||ii==57)
%         eiData(ii) = max(abs(temp.elecResp.cells.mainEI(ii,:))); %#ok<AGROW>
%     end
% end
clear temp



%% plotting pre-response probabilities as a function of stim pulse amplitude

plotColors(1,:) = [52 198 247]/255; %aqua

offsets(offsets == 0) = [];

%figure; hold on


for ii = 1:nPreElecs

    
    for kk = 1:length(prePulse(ii).amps)
        
        if prePulse(ii).analyzedPreResponse(kk) == 1
            
            figure; hold on
            
            for jj = 1:length(offsets)
                
                
                plotColor = (jj/length(offsets))*plotColors(1,:);
                
                
                stimAmps = abs(prePulse(ii).stimAmps{jj,kk});
                preRespRates = zeros(1, length(prePulse(ii).preResp{jj,kk}));
                for mm = 1:length(prePulse(ii).preResp{jj,kk})
                    preRespRates(mm) = sum(prePulse(ii).preResp{jj,kk}{mm})/length(prePulse(ii).preResp{jj,kk}{mm});
                end
                
                %fit a line to the data
                A = ones(length(stimAmps),2);
                A(:,1) = stimAmps;
                b = preRespRates';
                x = A\b;
                
                plot(stimAmps, x(1)*stimAmps + x(2), 'color', plotColor, 'lineWidth', 2)
                
                
            end
            hold off
            set(gca, 'ylim', [0 1])
        end
    end
end






