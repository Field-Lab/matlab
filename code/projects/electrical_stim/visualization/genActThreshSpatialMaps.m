function [allEiAmps, allVals, tempoutput] = genActThreshSpatialMaps(fPath,neuronId, varargin)
% GENACTTHRESHSPATIALMAPS(...) generates images that represent the
% activation thresholds of a cell for multiple stimulating electrodes a
% inputs:   fpath       path to elecResp files
%           neuronId    neuron id, can be a single value or a vector of
%                           multiple ids
%   optional:  saveFiles
%              printElecs
%              plotResponseCurves plots the sigmoidal activation curve for
%              each electrode that activeates a cell. could generate a lot
%              of figures
% outputs:  allEiAmps   matrix with ei info for all cells
%           allVals     matrix with activation info for all cells
% usage: genActThreshSpatialMaps('/Volumes/Analysis/2014-08-20-1/data003/', [7652 7473 7443 7427 7321 4366 4066 4007 3950 3946 1202 1096 752 544 256])
% genActThreshSpatialMaps('/Volumes/Analysis/2012-09-24-3/data004/', 4058,'printElecs',1);
%
% Lauren Grosberg 9/2014


% Set up defaults for optional parameters
printElecs = 0 ;
saveFiles = 1;
plotResponseCurves = 0;

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
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

dirInfo = dir(fPath);
for nn = 1:length(neuronId);
    clear eiamps elVals actThresholds
    neuron = neuronId(nn);
    for  n = 3:size(dirInfo,1) %%106 %
        if ~dirInfo(n).isdir
            fname = dirInfo(n).name;
            i = find(fname=='_',2,'first');
            if strcmp(['n' num2str(neuron)],fname(i(1)+1:i(2)-1))
                temp = load([fPath fname]);
                elecResp = temp.elecResp; clear temp;
                responseProb = elecResp.analysis.successRates;
                stimAmps     = abs(elecResp.stimInfo.stimAmps);
                
                % Define function that will be used to fit data
                % (F is a vector of fitting parameters)
                f = @(F,x) (1 +exp(-F(1)*(x - F(2)))).^(-1); % sigmoid
                F_fitted = nlinfit(stimAmps,responseProb,f,[1 1]);
                
                % Plot data fit
                y = f(F_fitted,stimAmps);
                
                % Find the threshold voltage for 50% response probability
                ii = find(y>0.5, 1,'first');
                if ii<size(stimAmps,1)
                    xx = stimAmps(ii-1):0.001:stimAmps(ii+1);
                    yy = f(F_fitted,xx);
                    
                    threshVoltage = xx(find(yy>0.5,1,'first'));
                    actThresholds(elecResp.stimInfo.electrodes) = threshVoltage;
                    
                    %Plotting
                    if plotResponseCurves
                        figure; plot(stimAmps, responseProb,'.'); ylim([0 1])
                        xlabel('Current (\muA)'); ylabel('response probability');
                        title(sprintf('neuron %d; stimulating electrode %d',neuron,elecResp.stimInfo.electrodes));
                        
                        hold on; plot(stimAmps,y,'g');
                        hold on; plot(threshVoltage,yy(find(yy>0.5,1,'first')),'ro');
                        text(threshVoltage+0.1,0.5,[num2str(threshVoltage) '\muA'])
                    end
                end
            end
        end
    end
    %
    elVals = zeros([1 512]);
    try
    elVals(1:size(actThresholds,2)) = actThresholds;
    catch
        keyboard; 
    end
    [~,aL_view] = ei2matrix(elVals);
    allVals(nn,:) = elVals;
    % figure; imagesc(1./aL_view,[0 1.2]); axis image; colorbar; colormap jet;
    %%
    figure;
    subplot(2,1,1); imagesc(aL_view); axis image; colorbar; colormap hot; axis off;
    % set(gcf,'Position',[25  300  1230 570]);
    title(['50% activation threshold n' num2str(neuron)]);
    tempoutput = elecResp.cells.mainEI; 
    eiamps = max(elecResp.cells.mainEI') - min(elecResp.cells.mainEI');
    [~,eipic] = ei2matrix(eiamps);
    subplot(2,1,2); imagesc(eipic,[0 max(eiamps)/2]); axis image; colorbar; colormap hot;
    title(['EI for n ' num2str(neuron)]); axis off;
    
    allEiAmps(nn,:) = eiamps;
    
    
    if printElecs
        h = load('/Users/grosberg/matlab/array_matrix_id510'); % Find a more general location for this or call a different text file.
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
end % Loop through neurons
end % end function