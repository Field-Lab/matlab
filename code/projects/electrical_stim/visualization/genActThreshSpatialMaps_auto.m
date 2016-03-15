function [allVals, actProb, actAmp] = genActThreshSpatialMaps_auto(fPath,eipath,neuronId, varargin)
% GENACTTHRESHSPATIALMAPS(...) generates images that represent the
% activation thresholds of a cell for multiple stimulating electrodes
% for use with automatic analysis
% inputs:   fpath       path to elecRespAuto files
%           neuronId    neuron id, can be a single value or a vector of
%                           multiple ids
%           eipath    path to ei file in order to load eis for plotting.
%           (must end in .ei)
%   optional: saveFiles
%             printElecs
%             plotResponseCurves: plots the sigmoidal activation curve for
%              each electrode that activeates a cell. could generate a lot
%              of figures
%             dispCheckedElecs: shows which electrodes were analyzed
%             assignValToPatternNo: useful in local return experiments, for
%              example, when threshold value is better assigned to a
%              primary electrode (here defined by the pattern number) than
%              to all the electrodes used in the stimulation pattern.
%              Default 0, set to 1
%             suppress plots: disables plotting; function returns outputs
%              only. Default 0.
%
% outputs:  allEiAmps   matrix with ei info for all cells
%           allVals     matrix with activation info for all cells
% usage: genActThreshSpatialMaps('/Volumes/Analysis/2014-08-20-1/data003/', [7652 7473 7443 7427 7321 4366 4066 4007 3950 3946 1202 1096 752 544 256])
% genActThreshSpatialMaps('/Volumes/Analysis/2012-09-24-3/data004/', 4058,'printElecs',1);
%
% Lauren Grosberg 12/2015

% Suppress warning: The Jacobian at the solution is ill-conditioned, and 
% some model parameters may not be estimated well (they are not identifi-
% able).  Use caution in making predictions.
warning('off', 'stats:nlinfit:IllConditionedJacobian');

% Set up defaults for optional parameters
printElecs = 0 ;
saveFiles = 0;
plotResponseCurves = 0;
dispCheckedElecs = 0; 
assignValToPatternNo = 0; 
suppressPlots = 0; 

nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn', ...
        'Unexpected number of arguments');
    throw(err);
end

% Reading the optional input arguments
for kk=1:(nbin/2)
    if ~ischar(varargin{kk*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{kk*2-1})
        case 'printelecs'
            printElecs = varargin{kk*2};
        case 'savefiles'
            saveFiles = varargin{kk*2};
        case 'plotresponsecurves'
            plotResponseCurves = varargin{kk*2};
        case 'dispcheckedelecs'
            dispCheckedElecs = varargin{kk*2}; 
        case 'assignvaltopatternno'
            assignValToPatternNo = varargin{kk*2};
        case 'suppressplots'
            suppressPlots = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

dirInfo = dir(fPath);
for nn = 1:length(neuronId);
    clear eiamps elVals % actThresholds
    neuron = neuronId(nn);
    for  n = 3:size(dirInfo,1) %%106 %
        if ~dirInfo(n).isdir && ~strcmp(dirInfo(n).name,'activationResults.mat')
            fname = dirInfo(n).name;
            i = find(fname=='_',2,'first');
            if strcmp(['n' num2str(neuron)],fname(i(1)+1:i(2)-1))
                temp = load([fPath fname]);
%                 disp([fPath fname]);
                elecResp = temp.elecRespAuto; clear temp;
                try
                    spikes = cell2mat(elecResp.spikes);
%                     spikes(find(isnan(spikes))) = 0; 
                    responseProb = sum(spikes,2)/size(spikes,2);
                    stimAmps     = abs(elecResp.stimInfo.listAmps);
                catch err
                    disp([fPath fname ' is empty']);
                    throw(err)
                end
                % LG temp(?) add 5/19/15 Find the real response probability >50%
                idx = find(responseProb>0.5,1,'first');
                if idx < size(stimAmps,1)
                    actProb(elecResp.stimInfo.patternNo) = responseProb(idx);
                    actAmp(elecResp.stimInfo.patternNo)  = stimAmps(idx);
                else
                    actProb(elecResp.stimInfo.patternNo) = 0;
                    actAmp(elecResp.stimInfo.patternNo) = max(stimAmps);
                end
                responseProb(isnan(responseProb)) =0; 
                [threshold, ~, ~] = fitToErf_inline(elecResp,...
                    responseProb,plotResponseCurves); 
                if threshold > 0
                    
                    if assignValToPatternNo
                        actThresholds(elecResp.stimInfo.patternNo) = threshold;
                    else
                        actThresholds(unique(elecResp.stimInfo.listStimElecs)) = threshold;
                    end
                else % Case where pattern did not activate the cell
                    if dispCheckedElecs
                        if assignValToPatternNo
                            actThresholds(elecResp.stimInfo.patternNo) = 10;
                        else
                            actThresholds(unique(elecResp.stimInfo.listStimElecs)) = 10;
                        end
                    end
                end
                
            end % end load elecRespAuto files from a particular neuron 
        end
    end % End search through directory
    %
    elVals = zeros([1 512]);
    try
        elVals(1:size(actThresholds,2)) = actThresholds;
    catch
        fprintf('no elecResp files for neuron %0.0f in %s\n',neuronId,fPath);
        allEiAmps = []; 
        allVals =zeros(1,512);  
        actProb = []; 
        actAmp = []; 
        return; 
    end
    clear actThresholds;
    [~,aL_view] = ei2matrix(elVals);
    allVals(nn,:) = elVals; %#ok<AGROW>
    [eiM,nIdList] = convertEiFileToMatrix(eipath);
    alleiamps = squeeze(max(eiM,[],2) - min(eiM,[],2)); 
    nIdx = find(neuron == nIdList);
    eiamps = alleiamps(:,nIdx);  %#ok<FNDSB>
    allEiAmps(nn,:) = eiamps'; %#ok<AGROW>
    if ~suppressPlots
        [xc,yc] = getElectrodeCoords512();
        figure;
        idx = find(eiamps>2);
        scatter(xc(idx),yc(idx),eiamps(idx),'k','filled');   hold on;
        [axon_x, axon_y, ~, soma_x, soma_y] = weighted_axon_poly_reg(eiamps);
        plot(axon_x, axon_y, 'Color', 'r','LineWidth',2, 'DisplayName', num2str(neuron));
        scatter(soma_x,soma_y,max(eiamps)+0.1,'r','filled');
        idx = find(elVals);
        if ~isempty(idx)
            scatter(xc(idx),yc(idx),100,elVals(idx),'filled');
        end
        axis image; axis off; c=colorbar; caxis([0.3 3]); 
        c.Label.String = 'activation threshold \muA';
        title(['EI for n ' num2str(neuron) ' -- 50% activation thresholds overlayed']); axis off;
       
        if printElecs
            for eIdx = 1:length(idx)
                text(xc(idx(eIdx)),yc(idx(eIdx))+10,num2str(idx(eIdx)),...
                    'HorizontalAlignment','center'); 
            end
        end
        
        if saveFiles
            [savename,savepath] = uiputfile('*.*','save images as');
            savingName = [savepath savename];
            saveas(gcf, savingName,'epsc');
            saveas(gcf, savingName,'jpg');
        end
    end %End if statement to show plots
    
end % Loop through neurons

    function [threshold,completeFit, erfErr] = fitToErf_inline(elecResp,responseProb,plotResponseCurves)
        nMovies = length(elecResp.stimInfo.listAmps);
        
        data = zeros(2, nMovies);
        data(1,:) = elecResp.stimInfo.listAmps;
        data(2,:) = responseProb;
        data(3,:) = elecResp.tracesInfo.I;
      
             
        % linear-based
        data(1,:) = abs(data(1,:));
        [erfParams, completeFit, erfErr] = erfFitter(data, 2, -1, 'makePlot', plotResponseCurves);
        threshold = -erfParams(2)/erfParams(1);
        
    end 
end% end function