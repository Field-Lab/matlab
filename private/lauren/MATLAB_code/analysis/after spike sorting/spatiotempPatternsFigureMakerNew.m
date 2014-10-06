function spatiotempPatternsFigureMakerNew(pathToData, neuronIDs, electrodes, patternOrders, offset, orderNo, movieNo, varargin)

% reworked and tested with spatiotemporal_2011-06-24-5_data010.m October
% 2011

% neuronIDs        neuron targets should correspond with pattern numbers (matches colors)
% electrodes       stim electrodes: order should correspond with target neurons in neuronIDs
% patternOrders    orders in which corresponding patterns are played (cell array index refers 
%                  to which sequence is applied, vector index corresponds to 'electrodes',
%                  value = ith pulse applied in sequence)
%                  * only the orderNoth order is plotted
% offset           offset between pulses in sequence, in ms (ignored if orderNo == 0 or 1)
% orderNo          which order in patternOrders to plot (cell array index)
%                  if orderNo == 0, plot simultaneous stim responses
%                  if orderNo == -1, plot isolated stim responses

p = inputParser;

p.addRequired('pathToData', @ischar)
p.addRequired('neuronIDs', @isnumeric)
p.addRequired('electrodes', @isnumeric)
p.addRequired('patternOrder', @iscell)
p.addRequired('offset', @isnumeric)
p.addRequired('orderNo', @isnumeric) %-1 signifies isolated stim, 0 signifies simultaneous stim, >0 signifies one of the sequences
p.addRequired('movieNo', @isnumeric)

colorSchemeDef(1,:) = [27 117 187]/256;
colorSchemeDef(2,:) = [190 30 45]/256;
colorSchemeDef(3,:) = [41 180 115]/256;
colorSchemeDef(4,:) = [0.8 0.5 0];

p.addParamValue('colorScheme', colorSchemeDef, @isnumeric)
p.addParamValue('xLimits', [], @isnumeric)
p.addParamValue('pxPerMs', 50, @isnumeric) %width of axes in pixels corresponding to 1 ms
p.addParamValue('rasterMarkerSize', 5, @isnumeric)

p.parse(pathToData, neuronIDs, electrodes, patternOrders, offset, orderNo, movieNo, varargin{:})

colorScheme = p.Results.colorScheme;
xLimits = p.Results.xLimits;
pxPerMs = p.Results.pxPerMs;
rasterMarkerSize = p.Results.rasterMarkerSize;

if orderNo > 0 || orderNo == -1 %not simultaneous stim
    nStim = length(electrodes);
    if orderNo > 0
        patternOrder = patternOrders{orderNo};
        if nStim ~= length(patternOrder)
            warning('this code is set to work only if each electrode corresponds to a single pulse in the pulse order. aborting!');
            return
        end
    end
end


if isempty(xLimits)
    if orderNo > 0 %not simultaneous or isolated stim
        xLimits = [0 offset*length(patternOrder)];
    else
        xLimits = [0 2.5];
    end
end

prevDir = pwd;

cd(pathToData)

warndlg(['note: timing of pulse must be checked -- sample at which pulse is assumed to be sample BEFORE current is flowing'...
    '(1 sample of amplifier disconnection before pulse)--need to check with Pawel whether this assumption is accurate'])

%% extract data

if orderNo > 0 || orderNo == -1 %not simultaneous stim
    latencies = cell(length(neuronIDs), nStim); %stores latency information from elecResps
    pulses = cell(nStim, 1);
    pulseAmps = zeros(nStim, 1);
    validReg = cell(length(neuronIDs), nStim);
    %electrodesAll = cell(nStim, 1);
    for jj = 1:nStim % jj = stim electrode
        for ii = 1:length(neuronIDs)
            if orderNo > 0 %sequential stim
                elecRespFile = ['elecResp_n' num2str(neuronIDs(ii)) '_p' num2str(electrodes(jj)) '_d' num2str(20*offset) '_o' num2str(orderNo) '.mat'];
            else %orderNo == -1, isolated stim
                elecRespFile = ['elecResp_n' num2str(neuronIDs(ii)) '_p' num2str(electrodes(jj)) '.mat'];
            end
            load(elecRespFile)
            movieInd = find(elecResp.stimInfo.movieNos == movieNo);

            if ~elecResp.analysis.finalized(movieInd)
                disp(['warning: movie ' num2str(movieNo) ' of ' elecRespFile ' is not locked...aborting!'])
                return
            end
            
            latencies{ii,jj} = elecResp.analysis.latencies{movieInd};
            validReg{ii,jj} = convertTempOffsetWindow(elecResp, movieInd);
        end
        pulses{jj} = elecResp.stimInfo.pulseVectors{movieInd}; % elec x time (1) vs. amp (2) x samples
        pulseAmps(jj) = elecResp.stimInfo.stimAmps(movieInd);
        %electrodesAll{jj} = elecResp.stimInfo.electrodes;
    end
    %stitch together latencies and valid regions
    if orderNo > 0
        latCat = cell(length(neuronIDs), 1);
        validRegCat = cell(length(neuronIDs), 1);
        for ii = 1:length(neuronIDs)
            latCat{ii} = cell(size(latencies{ii,1},1),1);
            for jj = 1:nStim
                for kk = 1:size(latencies{ii,jj},1) %trials
                    latTmp = latencies{ii,jj}(kk,:);
                    latCat{ii}{kk} = [latCat{ii}{kk} latTmp(latTmp~=0) + (patternOrder(jj)-1)*20*offset];
                end
                validRegCat{ii} = [validRegCat{ii}; (patternOrder(jj)-1)*20*offset + validReg{ii,jj}];
            end
        end
    end
else % orderNo == 0 %simultaneous stim
    latencies = cell(length(neuronIDs), 1);
    validReg = cell(length(neuronIDs), 1);
    for ii = 1:length(neuronIDs)
        elecRespFile = ['elecResp_n' num2str(neuronIDs(ii)) '_psim.mat'];
        load(elecRespFile)
        movieInd = find(elecResp.stimInfo.movieNos == movieNo);

        if ~elecResp.analysis.finalized(movieInd)
            disp(['warning: movie ' num2str(movieNo) ' of ' elecRespFile ' is not locked...aborting!'])
            return
        end
        latencies{ii} = elecResp.analysis.latencies{movieInd};

        validReg{ii} = convertTempOffsetWindow(elecResp, movieInd);
        pulses = elecResp.stimInfo.pulseVectors{movieInd};
    end
    electrodesSim = elecResp.stimInfo.electrodes;
end


%% shift latency and valid region values 1 sample left because first sample
% should be time = 0 (when stimulus is applied), but in elecResp files, samples start at time = 1 sample

for ii = 1:length(neuronIDs)
    if orderNo > 0 %sequential
        validRegCat{ii} = validRegCat{ii} - 1;
        for jj = 1:length(latCat{ii})
            latCat{ii}{jj} = latCat{ii}{jj} - 1;
        end
    elseif orderNo == -1 %isolated
        validReg{ii} = validReg{ii} - 1;
        for jj = 1:nStim
            validReg{ii,jj} = validReg{ii,jj} - 1;
            latencies{ii,jj}(latencies{ii,jj}~=0) = latencies{ii,jj}(latencies{ii,jj}~=0) - 1;
        end
    else %simultaneous
        validReg{ii} = validReg{ii} - 1;
        latencies{ii}(latencies{ii}~=0) = latencies{ii}(latencies{ii}~=0) - 1;
    end
end


%% check for refractory period violations (refractory period of 1 ms)

if orderNo > 0 %sequential
    for ii = 1:length(neuronIDs)
        for jj = 1:length(latCat{ii})
            latTmp = latCat{ii}{jj};
            latTmp = sort(latTmp);
            if any(diff(latTmp) < 0.5)
                warndlg('There is a refractory period (0.5 ms) violation in this data')
                keyboard
            end
        end
    end
end

%% make figure


if orderNo > 0 %sequential stim
    figure('position', [100 100 pxPerMs*xLimits(2)+50 length(neuronIDs)*60+100])
    
    %plot pulses
    pulseAxes = axes('units', 'pixels', 'position', [20 25 pxPerMs*xLimits(2) 60]);
    labelAxes = axes('units', 'pixels', 'position', [pxPerMs*xLimits(2)+20 25 30 60]);
    
    maxPulseRange = 0;
    for jj = 1:nStim
        maxPulseRange = max(maxPulseRange, max(pulses{jj}(1,2,:)) - min(pulses{jj}(1,2,:)));
    end
    
    for jj = 1:nStim
        axes(pulseAxes); hold on
        pulseVector = reshape(pulses{jj}(1,:,:), 2, []);
        plot([0, pulseVector(1,:)/1000 + (patternOrder(jj)-1)*offset, xLimits(2)],...
            [0, pulseVector(2,:)/(maxPulseRange*1.1), 0] + (nStim-jj)+1, 'Color', colorScheme(jj,:));
        hold off
        axes(labelAxes); hold on
        text(0.1, (nStim-jj)+1, num2str(electrodes(jj))); hold off
    end
    pulseYLims = [0 nStim + 1];

    set(pulseAxes, 'xLim', xLimits, 'yLim', pulseYLims)
    set(labelAxes, 'yLim', get(pulseAxes, 'yLim'))
    axes(pulseAxes); xlabel('time (ms)')
    axes(labelAxes); axis off
    
    %plot rasters
    for ii = 1:length(neuronIDs)
        axes('units', 'pixels', 'position', [20 (length(neuronIDs) - ii)*60 + 100 pxPerMs*xLimits(2) 50]); %entire row
        hold on
                
        %grey out unanalyzed regions
        valRegTmp = sort(validRegCat{ii}/20);
        yReg = [0 0 length(latCat{ii})+1 length(latCat{ii})+1];
        patch([0 valRegTmp(1,1) valRegTmp(1,1) 0], yReg, [0.9 0.9 0.9], 'edgeColor', 'none')
        for kk = 2:size(valRegTmp,1)
            patch([valRegTmp(kk-1,2) valRegTmp(kk,1) valRegTmp(kk,1) valRegTmp(kk-1,2)], yReg, [0.9 0.9 0.9], 'edgeColor', 'none')
        end
        patch([valRegTmp(end,2) xLimits(2) xLimits(2) valRegTmp(end,2)], yReg, [0.9 0.9 0.9], 'edgeColor', 'none')

        %plot raster points
        for kk = 1:length(latCat{ii}) %trials
            plot(latCat{ii}{kk}/20, kk*ones(size(latCat{ii}{kk})), 'o', 'markerSize', rasterMarkerSize, 'markerFaceColor', colorScheme(ii,:), 'markerEdgeColor', 'none')
        end
        
        %plot box
        plot([0 xLimits(2) xLimits(2) 0 0], [0 0 length(latCat{ii})+1 length(latCat{ii})+1 0], 'k-')
        hold off

        
        set(gca, 'xLim', xLimits, 'yLim', [0 length(latCat{ii})+1],'box', 'on',...
            'xTickLabel', [], 'yTickLabel', [])
        
        %plots neuron IDs
        axes('units', 'pixels', 'position', [pxPerMs*xLimits(2)+20 (length(neuronIDs) - ii)*60+100 30 50])
        axis off
        text(0.1, 0.4, num2str(neuronIDs(ii)))
        set(gca, 'xlim', [0 1])
    end
    

elseif orderNo == 0 %simultaneous stim
    figure('position', [100 100 pxPerMs*xLimits(2)+50 length(neuronIDs)*60+100])
    
    %plot pulses
    pulseAxes = axes('units', 'pixels', 'position', [20 25 pxPerMs*xLimits(2) 60]);
    labelAxes = axes('units', 'pixels', 'position', [pxPerMs*xLimits(2)+20 25 30 60]);
    
    maxPulseRange = max(max(pulses(:,2,:))) - min(min(pulses(:,2,:)));
    
    for jj = 1:length(electrodes)
        axes(pulseAxes); hold on
        pulseVector = reshape(pulses(electrodesSim == electrodes(jj),:,:), 2, []);
        plot([0, pulseVector(1,:)/1000, xLimits(2)],...
            [0, pulseVector(2,:)/(maxPulseRange*1.1), 0] + (length(electrodes)-jj)+1, 'Color', colorScheme(jj,:));
        hold off
        axes(labelAxes); hold on
        text(0.1, (length(electrodes)-jj)+1, num2str(electrodes(jj))); hold off
    end
    
    pulseYLims = [0 length(electrodes) + 1];

    set(pulseAxes, 'xLim', xLimits, 'yLim', pulseYLims)
    set(labelAxes, 'yLim', get(pulseAxes, 'yLim'))
    axes(pulseAxes); xlabel('time (ms)')
    axes(labelAxes); axis off
    
    %plot rasters
    for ii = 1:length(neuronIDs)
        axes('units', 'pixels', 'position', [20 (length(neuronIDs) - ii)*60 + 100 pxPerMs*xLimits(2) 50]); %entire row
        hold on
        
        %grey out unanalyzed regions
        yReg = [0 0 length(latencies{ii})+1 length(latencies{ii})+1];
        patch([0 validReg{ii}(1)/20 validReg{ii}(1)/20 0], yReg, [0.9 0.9 0.9], 'edgeColor', 'none')
        patch([validReg{ii}(2)/20 xLimits(2) xLimits(2) validReg{ii}(2)/20], yReg, [0.9 0.9 0.9], 'edgeColor', 'none')

        
        %plot raster points
        for kk = 1:length(latencies{ii}) %trials
            if latencies{ii}(kk) ~= 0
                plot(latencies{ii}(kk)/20, kk, 'o', 'markerSize', rasterMarkerSize, 'markerFaceColor', colorScheme(ii,:), 'markerEdgeColor', 'none')
            end
        end
        
        %plot box
        plot([0 xLimits(2) xLimits(2) 0 0], [0 0 length(latencies{ii})+1 length(latencies{ii})+1 0], 'k-')
        hold off
        
        
        set(gca, 'xLim', xLimits, 'yLim', [0 length(latencies{ii})+1],'box', 'on',...
            'xTickLabel', [], 'yTickLabel', [])

        %plots neuron IDs
        axes('units', 'pixels', 'position', [pxPerMs*xLimits(2)+20 (length(neuronIDs) - ii)*60+100 30 50])
        axis off
        text(0.1, 0.4, num2str(neuronIDs(ii)))
        set(gca, 'xlim', [0 1])
    end
    
else %orderNo == -1 %isolated stim
    figure('position', [100 100 pxPerMs*xLimits(2)*nStim + 20*nStim + 50 length(neuronIDs)*60+100])
    
    
    %plot pulses
    labelAxes = axes('units', 'pixels', 'position', [nStim*xLimits(2)*pxPerMs + nStim*20 25 30 60]);
    pulseYLims = [0 nStim + 1];

    maxPulseRange = 0;
    for jj = 1:nStim
        maxPulseRange = max(maxPulseRange, max(pulses{jj}(1,2,:)) - min(pulses{jj}(1,2,:)));
    end
    
    pulseAxes = cell(nStim,1);
    for jj = 1:nStim
        pulseAxes{jj} = axes('units', 'pixels', 'position', [20+(jj-1)*xLimits(2)*pxPerMs + (jj-1)*20 25 pxPerMs*xLimits(2) 60]);
        hold on
        for kk = 1:length(electrodes)
            if jj == kk
                pulseVector = reshape(pulses{jj}(1,:,:), 2, []);
            else
                pulseVector = [0; 0];
            end
            plot([0, pulseVector(1,:)/1000, xLimits(2)],...
                [0, pulseVector(2,:)/(maxPulseRange*1.1), 0] + (nStim-kk)+1, 'Color', colorScheme(kk,:));
            set(pulseAxes{jj}, 'xLim', xLimits, 'yLim', pulseYLims)
            xlabel('time (ms)')
        end
        hold off
        axes(labelAxes); hold on
        text(0.1, (nStim-jj)+1, num2str(electrodes(jj))); hold off
    end    
    set(labelAxes, 'yLim', get(pulseAxes{1}, 'yLim'))
    axes(labelAxes); axis off
    
    %plot rasters
    for ii = 1:length(neuronIDs)
        for jj = 1:nStim
            axes('units', 'pixels', 'position',...
                [20+(jj-1)*xLimits(2)*pxPerMs+(jj-1)*20 (length(neuronIDs)-ii)*60+100 pxPerMs*xLimits(2) 50]);
            hold on
            
            %grey out unanalyzed regions
            yReg = [0 0 length(latencies{ii,jj})+1 length(latencies{ii,jj})+1];
            patch([0 validReg{ii,jj}(1)/20 validReg{ii,jj}(1)/20 0], yReg, [0.9 0.9 0.9], 'edgeColor', 'none')
            patch([validReg{ii,jj}(2)/20 xLimits(2) xLimits(2) validReg{ii,jj}(2)/20], yReg, [0.9 0.9 0.9], 'edgeColor', 'none')
            
            
            %plot raster points
            for kk = 1:length(latencies{ii,jj}) %trials
                if latencies{ii,jj}(kk) ~= 0
                    plot(latencies{ii,jj}(kk)/20, kk, 'o', 'markerSize', rasterMarkerSize, 'markerFaceColor', colorScheme(ii,:), 'markerEdgeColor', 'none')
                end
            end
            
            %plot box
            plot([0 xLimits(2) xLimits(2) 0 0], [0 0 length(latencies{ii,jj})+1 length(latencies{ii,jj})+1 0], 'k-')
            hold off
            
            set(gca, 'xLim', xLimits, 'yLim', [0 length(latencies{ii,jj})+1],'box', 'on',...
                'xTickLabel', [], 'yTickLabel', [])
        end
                    
        %plots neuron IDs
        axes('units', 'pixels', 'position', [20+nStim*xLimits(2)*pxPerMs+(nStim-1)*20 (length(neuronIDs) - ii)*60+100 30 50])
        axis off
        text(0.1, 0.4, num2str(neuronIDs(ii)))
        set(gca, 'xlim', [0 1])
        
    end
end


cd(prevDir)

