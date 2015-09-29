function [allEiAmps, allVals, actAmp, actThreshold] = genActThreshSpatialMaps_Gonzalo(fPath,neuronId, Outputs,varargin)
% GENACTTHRESHSPATIALMAPS(...) generates images that represent the
% activation thresholds of a cell for multiple stimulating electrodes a
% inputs:   fpath       path to elecResp files
%           neuronId    neuron id, can be a single value or a vector of
%                           multiple ids
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
%
% outputs:  allEiAmps   matrix with ei info for all cells
%           allVals     matrix with activation info for all cells
% usage: genActThreshSpatialMaps('/Volumes/Analysis/2014-08-20-1/data003/', [7652 7473 7443 7427 7321 4366 4066 4007 3950 3946 1202 1096 752 544 256])
% genActThreshSpatialMaps('/Volumes/Analysis/2012-09-24-3/data004/', 4058,'printElecs',1);
%
% Lauren Grosberg 9/2014

% Suppress warning: The Jacobian at the solution is ill-conditioned, and
% some model parameters may not be estimated well (they are not identifi-
% able).  Use caution in making predictions.
warning('off', 'stats:nlinfit:IllConditionedJacobian');

% Set up defaults for optional parameters
printElecs = 0 ;
saveFiles = 0;
plotResponseCurves = 0;
dispCheckedElecs = 0;
assignValToPatternNo = 1;
suppressPlots = 0;
actAmp = zeros(2,512);
actThreshold = zeros(2,512);

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
    
    
    for j=1:length(Outputs)
        neuronO=-1;
        try neuronO=Outputs(j).neuronInfo.neuronIds;
        end
        
        if(isequal(neuronO,neuron))
            Output=Outputs(j);
            pathToAnalysisData = Output.path;
            patternNo = Output.stimInfo.patternNo;
            neuronIds = Output.neuronInfo.neuronIds;
            tempPath = [pathToAnalysisData 'elecResp_n' num2str(neuronIds) '_p' num2str(patternNo) '.mat'];
            temp     = load(tempPath);
            elecResp = temp.elecResp;
            
            [thresholdHum thresholdAlg] = fitToErfOutputAndHuman(Output);
            try
                responseProb = elecResp.analysis.successRates;
                stimAmps     = abs(elecResp.stimInfo.stimAmps);
            catch err
                disp([fPath fname ' is empty']);
                throw(err)
            end
            
            actAmp(1:2,elecResp.stimInfo.patternNo)=[thresholdHum thresholdAlg]';
            
            % L temp(?) add 5/19/15 Find the real response probability >50%
            %idx = find(responseProb>0.5,1,'first');
            
            %if idx < size(stimAmps,1)
             %   actProb(elecResp.stimInfo.patternNo) = responseProb(idx);
              %  actAmp(1:2,elecResp.stimInfo.patternNo)  = stimAmps(idx);
            %else
             %   actProb(elecResp.stimInfo.patternNo) = 0;
              %  actAmp(elecResp.stimInfo.patternNo) = max(stimAmps);
           % end
            
            if 1 > 0
                
                if assignValToPatternNo
                    actThresholds(1:2,elecResp.stimInfo.patternNo) =[thresholdHum thresholdAlg]';
               
                else
                    actThresholds(1:2,elecResp.stimInfo.electrodes) = [thresholdHum thresholdAlg]';
                end
            else % Case where pattern did not activate the cell
                if dispCheckedElecs
                    if assignValToPatternNo
                        actThresholds(1:2,elecResp.stimInfo.patternNo) = [thresholdHum thresholdAlg]';
                    else
                        actThresholds(1:2,elecResp.stimInfo.electrodes) = [thresholdHum thresholdAlg]';
                    end
                end
            end
            
            %                 % Define function that will be used to fit data
            %                 % (F is a vector of fitting parameters)
            %                 f = @(F,x) (1 +exp(-F(1)*(x - F(2)))).^(-1); % sigmoid
            %                 F_fitted = nlinfit(stimAmps,responseProb,f,[1 1]);
            %
            %                 % Plot data fit
            %                 y = f(F_fitted,stimAmps);
            %
            %                 % Find the threshold voltage for 50% response probability
            %                 ii = find(y>0.5, 1,'first');
            %                 if ii<size(stimAmps,1)
            %                     try
            %                     xx = stimAmps(ii-1):0.001:stimAmps(ii+1);
            %                     catch
            %                         keyboard
            %                     end
            %                     yy = f(F_fitted,xx);
            %
            %                     threshVoltage = xx(find(yy>0.5,1,'first'));
            %                     if assignValToPatternNo
            %                         actThresholds(elecResp.stimInfo.patternNo) = threshVoltage;
            %                     else
            %                         actThresholds(elecResp.stimInfo.electrodes) = threshVoltage;
            %                     end
            %
            %                     %Plotting
            %                     if plotResponseCurves
            %                         figure; plot(stimAmps, responseProb,'.'); ylim([0 1])
            %                         xlabel('Current (\muA)'); ylabel('response probability');
            %                         title(sprintf('neuron %d; stimulating electrode %d',neuron,elecResp.stimInfo.electrodes));
            %
            %                         hold on; plot(stimAmps,y,'g');
            %                         hold on; plot(threshVoltage,yy(find(yy>0.5,1,'first')),'ro');
            %                         text(threshVoltage+0.1,0.5,[num2str(threshVoltage) '\muA'])
            %                     end
            %                 else
            %                     if dispCheckedElecs
            %                         actThresholds(elecResp.stimInfo.patternNo) = 10;
            %                     end
            %                 end
        % end load elecResp files from a particular neuron
    end
end
% End search through directory
%
    elVals = zeros([1 512]);
    try
        elVals(1:size(actThresholds,2)) = 1;
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
    eiamps = max(elecResp.cells.mainEI,[],2) - min(elecResp.cells.mainEI,[],2);

    allEiAmps(nn,:) = eiamps'; %#ok<AGROW>
    if ~suppressPlots
        figure;
        subplot(2,1,1); imagesc(aL_view); axis image; colorbar; colormap hot; axis off;
        % set(gcf,'Position',[25  300  1230 570]);
        title(['50% activation threshold n' num2str(neuron)]);
        [~,eipic] = ei2matrix(eiamps);
        subplot(2,1,2); imagesc(eipic,[0 max(eiamps)/2]); axis image; colorbar; colormap hot;
        title(['EI for n ' num2str(neuron)]); axis off;
       
        if printElecs
            fname = fullfile(fileparts(mfilename('fullpath')),'../resources/array_matrix_id510');
            h = load(fname);
            electrodeMatrix = h.array_matrix_id510; clear h;
            for x = 1:size(electrodeMatrix,2)
                for y = 1:size(electrodeMatrix,1)
                    if y/2 == round(y/2)
                        th=text(x*2-1,y*2,num2str(electrodeMatrix(y,x)));
                        set(th,'color','g');
                    else
                        th=text(x*2,y*2,num2str(electrodeMatrix(y,x)));
                        set(th,'color','g');
                    end
                end
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
end % end function