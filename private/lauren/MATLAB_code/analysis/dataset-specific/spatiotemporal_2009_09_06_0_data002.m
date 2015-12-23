clear all

%cd /snle/lab/Experiments/Array/Analysis/2009-09-06-0/data002
cd('/Volumes/Palace/Analysis/Lauren/2009-09-06-0/new elecResp')

neuronIDs  = [93 407 519 634;   93  407 519 634;  407 519 634 93;   634 519 407 93;   519 407 93  634];
patternNos = [1  1   1   1;     201 301 401 501;  302 402 502 202;  503 403 303 203;  404 304 204 504];
movieNo = 8;
nPatterns = length(patternNos);

patternTimes = [0 0 0 0;  0 1 2 3;  0 1 2 3;  0 1 2 3;  0 1 2 3];

colorScheme(1,:) = [27 117 187]/256;
colorScheme(2,:) = [190 30 45]/256;
colorScheme(3,:) = [41 180 115]/256;
colorScheme(4,:) = [0.8 0.5 0];


% latencies = cell(size(neuronIDs));
% pulseVectors = cell(size(neuronIDs));
% for i = 1:nPatterns
%     temp = load(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(i))]);
%     elecResp = temp.elecResp;
%     if ~elecResp.analysis.finalized(movieNo)
%         disp(['warning: movie ' num2str(movieNo) ' of pattern ' num2str(patternNos(i)) ' is not locked'))
%     end
%     movieIndex = find(elecResp.stimInfo.movieNos == movieNo);
%     
%     latenciesAllNeurons{i} = latencies{movieIndex};
%     
%     pulses{i} = squeeze(elecResp.stimInfo.pulseVectors{movieIndex}); % elec x time (1) vs. amp (2) x samples
% end


%% plotting

latencies = cell(size(neuronIDs));
pulses = cell(size(neuronIDs));
for j = 1:size(neuronIDs, 2) %subsequent pulses
    for i = 1:size(neuronIDs, 1) % subsequent pulse orders
        temp = load(['elecResp_n' num2str(neuronIDs(i,j)) '_p' num2str(patternNos(i,j))]);
        elecResp = temp.elecResp;
        if ~elecResp.analysis.finalized(movieNo)
            disp(['warning: movie ' num2str(movieNo) ' of pattern ' num2str(patternNos(i,j)) ' is not locked'])
        end
        movieIndex = find(elecResp.stimInfo.movieNos == movieNo);

        latencies{i,j} = elecResp.analysis.latencies{movieIndex};
        pulses{i,j} = elecResp.stimInfo.pulseVectors{movieIndex}; % elec x time (1) vs. amp (2) x samples
    end
end

save('plotInfo.mat')


for i = 1:size(neuronIDs, 1)
    figure('position', [100 100 400 500])
    
    % plotting pulses
    axes('position', [0.05 0.05 0.9 0.2])
    hold on
    for j = 1:size(neuronIDs, 2)
        elecID = round(patternNos(i,j)/100) - 1;
        if size(pulses{i,j}, 1) == 1 %single electrode
            pulseVector = reshape(pulses{i,j}(1, :, :), 2, []);
            current = plot([0 pulseVector(1,:)/1000+patternTimes(i,j) 4],[1.5*elecID pulseVector(2,:) + 1.5*elecID 1.5*elecID]);
            set(findobj(current,'Type','line'),'Color', colorScheme(elecID,:), 'LineWidth', 1)
        else %all electrodes
            pulseVector = reshape(pulses{i,j}(j, :, :), 2, []);
            current = plot([0 pulseVector(1,:)/1000+patternTimes(i,j) 4],[1.5*j pulseVector(2,:) + 1.5*j 1.5*j]);
            set(findobj(current,'Type','line'),'Color', colorScheme(j,:), 'LineWidth', 1)
        end
    end
    hold off
    set(gca, 'xLim', [0 4])
    xlabel('time (ms)')
    
    % plotting spikes
    for j = 1:size(neuronIDs, 2)
        if size(pulses{i,j}, 1) == 1 %single electrode
            elecID = round(patternNos(i,j)/100) - 1;
            axes('position', [0.05 0.275+(elecID-1)*.175 0.9 0.15])
            hold on
            for k = 1:length(latencies{i,j})
                if latencies{i,j}(k) ~= 0
                    plot(latencies{i,j}(k)/20+patternTimes(i,j), k, 'k.', 'markerSize', 5,...
                        'markerFaceColor', colorScheme(elecID,:), 'markerEdgeColor', colorScheme(elecID,:))
                end
            end
            plot([patternTimes(i,j) patternTimes(i,j)], [0 1])
            hold off
            set(gca, 'xLim', [0 4], 'yLim', [0 length(latencies{i,j})+1])
            axis off
        else
            elecID = round(patternNos(i,j)/100) - 1;
            axes('position', [0.05 0.275+(j-1)*0.175 0.9 0.15])
            hold on
            for k = 1:length(latencies{i,j})
                if latencies{i,j}(k) ~= 0
                    plot(latencies{i,j}(k)/20+patternTimes(i,j), k, 'k.', 'markerSize', 5,...
                        'markerFaceColor', colorScheme(j,:), 'markerEdgeColor', colorScheme(j,:))
                end
            end
            hold off
            set(gca, 'xLim', [0 4], 'yLim', [1 length(latencies{i,j})])
            axis off
        end
    end
end

%% hacked version of the third spatiotemporal pattern (patterns *03) that plots raster over entire
% window of time for each cell

neuronIDs = [93 407 519 634];
patternNos = [503 403 303 203];

%neuronIDs  = [93  407 519 634;  407 519 634 93;   634 519 407 93;   519 407 93  634];
%patternNos = [503 403 303 203;  503 403 303 203;  503 403 303 203;  503 403 303 203];
movieNo = 8;
nPatterns = length(patternNos);

%patternTimes = [0 1 2 3;  0 1 2 3;  0 1 2 3;  0 1 2 3];



latencies = cell(size(neuronIDs));
pulses = cell(length(patternNos), 1);
for j = 1:length(patternNos) % subsequent patterns
    for i = 1:length(neuronIDs) %subsequent neurons

        temp = load(['elecResp_n' num2str(neuronIDs(i)) '_p' num2str(patternNos(j))]);
        elecResp = temp.elecResp;
        if ~elecResp.analysis.finalized(movieNo)
            disp(['warning: movie ' num2str(movieNo) ' of pattern ' num2str(patternNos(i,j)) ' is not locked'])
        end
        movieIndex = find(elecResp.stimInfo.movieNos == movieNo);

        latencies{i,j} = elecResp.analysis.latencies{movieIndex};
    end
    pulses{j} = elecResp.stimInfo.pulseVectors{movieIndex}; % elec x time (1) vs. amp (2) x samples
end



figure('position', [100 100 400 500])
for i = 1:length(neuronIDs) %progresses through neurons
    
    % broken!!!!!
%     % plotting pulses
%     axes('position', [0.05 0.05 0.9 0.2])
%     hold on
%     for j = 1:length(patternNos) %different patterns
%         elecID = 5-j; %since patterns are listed in reverse order
%         pulseVector = reshape(pulses{j}(1, :, :), 2, []);
%         plot([0 pulseVector(1,:)/1000 + elecID-1 4], [1.5*elecID pulseVector(2,:) + 1.5*elecID 1.5*elecID],...
%             'Color', colorScheme(elecID,:));
%     end
%     hold off
%     set(gca, 'xLim', [0 4])
%     xlabel('time (ms)')

    % plotting spikes
    for j = 1:length(patternNos)
        axes('position', [0.05 0.275+(i-1)*.175 0.9 0.15]) %entire row
        hold on
        for k = 1:length(latencies{i,j})
            if latencies{i,j}(k) ~= 0
                plot(latencies{i,j}(k)/20+ j-1, k, '.', 'markerSize', 5,...
                    'markerFaceColor', colorScheme(i,:), 'markerEdgeColor', colorScheme(i,:))
            end
        end
        %plot([patternTimes(i,j) patternTimes(i,j)], [0 1])
        hold off
        set(gca, 'xLim', [0 4], 'yLim', [0 length(latencies{i,j})+1])
        axis off
    end
end





