function cellIDs = showCellsNearElec(datarun,electrodeNo,varargin)
% Function that plots RFs of cells with EIs passing over a stimulating
% electrode
% inputs: 
%       datarun (need load_data, load_params, load_sta, load_sta_fits,
%                   get_sta_fits_from_vision, load_ei)
%       electrodeNo electrode number over which you are looking for signal
%  (optional)
%       dataPath - path to elecResp files for plotting stimulating
%                   threshold data with a slider plot
%       autoFile - set to true to look for elecRespAuto files
%       threshold - threshold to include EIs (default is 7)
% Lauren Grosberg, modified by Alena Rott summer 2015
% LG edits added 2/2016
% Usage example: 
%   dataPath = '/Volumes/Analysis/2015-11-09-3/data000/data000'; 
%   autoFilePath = '/Volumes/Analysis/2015-11-09-3/data001-data002-autosort/';
%   datarun = load_data(dataPath);
%   datarun = load_neurons(datarun);
%   datarun = load_sta(datarun, 'load_sta', 'all');
%   datarun = load_params(datarun);
%   datarun = load_ei(datarun, 'all');
%   showCellsNearElec(datarun,200,'datapath',autoFilePath,'autoFile',true,'threshold',10);

if length(varargin) ~= 0
    dataPath = varargin{1};    
end
% Set up defaults for optional parameters
dataPath = [] ;
autoFile = false;
threshold = 7;

nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

% Read the optional input arguments
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'datapath'
            dataPath = varargin{j*2};
        case 'autofile'
            autoFile = varargin{j*2};
        case 'threshold'
            threshold = varargin{j*2};      
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

positions  = datarun.ei.position;
eiFile   = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(datarun.names.rrs_ei_path);
try [cellSpec, cellType] = get_cell_indices( datarun, {'ON parasol', 'OFF parasol', 'ON midget', 'OFF midget'} );
    
    % ON Parasol cells
    ONPcellIDs = [];
    cellTypeIndex = strcmp(cellType,'ON parasol');
    neurons  = datarun.cell_ids(cellSpec{cellTypeIndex}); % Analyze a subset of neurons
    for i = 1:1:size(neurons,2) % For all ON parasol cells
        tempID  = neurons(i);          % gets the id number of the ith neuron
        ei = eiFile.getImage(tempID);  % gets ei data for ith neuron, storing the information as a 3D array
        ei = squeeze(ei(1,2:end,:))';
        eiAmps = max(ei)-min(ei);
        if eiAmps(electrodeNo) > threshold
            ONPcellIDs = cat(1,ONPcellIDs,tempID);
        end
    end
    
    % OFF Parasol cells
    OFFPcellIDs = [];
    cellTypeIndex = strcmp(cellType,'OFF parasol');
    neurons  = datarun.cell_ids(cellSpec{cellTypeIndex}); % Analyze a subset of neurons
    for i = 1:1:size(neurons,2) % For all ON parasol cells
        tempID  = neurons(i);          % gets the id number of the ith neuron
        ei = eiFile.getImage(tempID);  % gets ei data for ith neuron, storing the information as a 3D array
        ei = squeeze(ei(1,2:end,:))';
        eiAmps = max(ei)-min(ei);
        if eiAmps(electrodeNo) > threshold
            OFFPcellIDs = cat(1,OFFPcellIDs,tempID);
        end
    end
    
    % ON midget cells
    ONMcellIDs = [];
    cellTypeIndex = strcmp(cellType,'ON midget');
    neurons  = datarun.cell_ids(cellSpec{cellTypeIndex}); % Analyze a subset of neurons
    for i = 1:1:size(neurons,2) % For all ON parasol cells
        tempID  = neurons(i);          % gets the id number of the ith neuron
        ei = eiFile.getImage(tempID);  % gets ei data for ith neuron, storing the information as a 3D array
        ei = squeeze(ei(1,2:end,:))';
        eiAmps = max(ei)-min(ei);
        if eiAmps(electrodeNo) > threshold
            ONMcellIDs = cat(1,ONMcellIDs,tempID);
        end
    end
    
    % OFF midget cells
    OFFMcellIDs = [];
    cellTypeIndex = strcmp(cellType,'OFF midget');
    neurons  = datarun.cell_ids(cellSpec{cellTypeIndex}); % Analyze a subset of neurons
    for i = 1:1:size(neurons,2) % For all ON parasol cells
        tempID  = neurons(i);          % gets the id number of the ith neuron
        ei = eiFile.getImage(tempID);  % gets ei data for ith neuron, storing the information as a 3D array
        ei = squeeze(ei(1,2:end,:))';
        eiAmps = max(ei)-min(ei);
        if eiAmps(electrodeNo) > threshold
            OFFMcellIDs = cat(1,OFFMcellIDs,tempID);
        end
    end
    cellIDs = sort(cat(1,ONPcellIDs, OFFPcellIDs, ONMcellIDs, OFFMcellIDs));
    
    % Plot overlay of RF fits and cell IDs
    fh=figure; set(gcf,'Color',[1 1 1],'Position',[40 545 1785 555]); subplot(1,2,1); hold on;
    plot_rf_summaries(datarun, ONPcellIDs, 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b');
    plot_rf_summaries(datarun, OFFPcellIDs, 'clear', false, 'label', true, 'label_color', 'k', 'plot_fits', true, 'fit_color', 'k');
    plot_rf_summaries(datarun, ONMcellIDs, 'clear', false, 'label', true, 'label_color', 'm', 'plot_fits', true, 'fit_color', 'm');
    plot_rf_summaries(datarun, OFFMcellIDs, 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'g');
    axis image; axis off;
    title(sprintf('ON parasol (blue) OFF parasol (black)\n ON midget (pink) OFF midget (green)\n stimulating electrode %0.0f',electrodeNo));
   
catch % Plot all cells, unclassified by type. 
    [cellSpec, cellType] = get_cell_indices( datarun, 'all' );
    cellIDs = [];
    
    neurons  = datarun.cell_ids(cellSpec); % Analyze neurons
    for i = 1:1:size(neurons,2) % For all ON parasol cells
        tempID  = neurons(i);          % gets the id number of the ith neuron
        ei = eiFile.getImage(tempID);  % gets ei data for ith neuron, storing the information as a 3D array
        ei = squeeze(ei(1,2:end,:))';
        eiAmps = max(ei)-min(ei);
        if eiAmps(electrodeNo) > threshold
            cellIDs = cat(1,cellIDs,tempID);
        end
    end
    % Plot RF fits and cell IDs
    fh=figure; set(gcf,'Color',[1 1 1],'Position',[40 545 1785 555]); subplot(1,2,1); hold on;
    plot_rf_summaries(datarun, cellIDs, 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b');
    axis image; axis off;
    % legend('ON parasol','OFF parasol','ON midget','OFF midget')
    title(sprintf('%s stimulating electrode %0.0f',datarun.names.short_name,electrodeNo));
    
end %end catch statement
eisToPlot = get_cell_indices(datarun,cellIDs);
neurons = datarun.cell_ids(:)';
%{
% Scatterplot showing EIs overlayed
figure; set(gcf,'Color',[1 1 1]); hold on;


maxRecElec = spalloc(size(neurons,2),1,length(eisToPlot));
for i = 1:length(eisToPlot)
    tempID  = neurons(eisToPlot(i));          % gets the id number of the ith neuron
    ei = eiFile.getImage(tempID);  % gets ei data for ith neuron, storing the information as a 3D array
    ei = squeeze(ei(1,2:end,:))';
    eiAmps = max(ei)-min(ei);
    eiAmps(find(eiAmps<threshold)) = 0;
    maxVal = max(eiAmps);
    maxRecElec(tempID) = find(eiAmps == maxVal,1);
    scatter(positions(:,1),positions(:,2),10*eiAmps+0.1,'filled');
end
axis image; axis off; 
text(positions(maxRecElec(cellIDs),1),positions(maxRecElec(cellIDs),2),...
    num2str(cellIDs),'HorizontalAlignment','center');
% Plot stimulating electrode
scatter(positions(electrodeNo,1),positions(electrodeNo,2),900,'k'); 
scatter(positions(electrodeNo,1),positions(electrodeNo,2),300,'k', 'filled');
% text(positions(electrodeNo,1),positions(electrodeNo,2),'stim elec');
title(sprintf('EIs with signal above %0.0f on electrode %0.0f',threshold,electrodeNo)); 
%}

figure(fh); subplot(1,2,2); hold on;
colors = hsv(length(eisToPlot));
e1 = scatter(positions(electrodeNo,1),positions(electrodeNo,2),900,'k'); 
e2 = scatter(positions(electrodeNo,1),positions(electrodeNo,2),300,'k', 'filled'); 
set(get(get(e1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(e2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

dist_to_elec = [];
neuronIDs = [];
%Can be easily edited to access data from an elecRespAuto file.
if ~isempty(dataPath)
    files = dir(fullfile(dataPath));
    k = 1;
    for n = 1:length(files)
        if ~files(n).isdir
            fileName = files(n).name;
            autoFile = strcmp(fileName(9:12),'Auto');
            % Only load pattern of interest by reading the filename. This
            % will have to be updated when patternNo ~= stimElectrode
            ii = find(fileName=='_');
            searchPattern = ['p' num2str(electrodeNo) '.'];
            if strcmp(searchPattern,fileName(ii(2)+1:ii(2)+length(searchPattern)))
                load(strcat(dataPath, fileName));
                if autoFile
                    neuronID = elecRespAuto.neuronInfo.neuronIds;
                else
                    neuronID = elecResp.cells.main;
                end
                if ismember(neuronID, cellIDs)
                    neuronIDs(k) = neuronID; % Adds many IDs
                    if autoFile
                        stim_levels{k, :} = elecRespAuto.stimInfo.listAmps(:,find(min(elecRespAuto.stimInfo.listAmps,[],1) == min(elecRespAuto.stimInfo.listAmps(:)))); %#ok<FNDSB>
                        activation{k, :} = elecRespAuto.LogisticReg;                    
                    else
                        stim_levels{k, :} = elecResp.stimInfo.stimAmps;
                        activation{k, :} = elecResp.analysis.successRates;
                    end
                    disp(k)
                    k = k + 1;
                end % End if statement requiring that the neuron is a member of the list of cellIDs
            end % End if statement requiring that the patternNo is the same as the stimulating electrode
        end % End requirement that elecResp file is not a directory
    end % End search through files in elecResp(Auto) directory
end

for i = 1:length(eisToPlot) %number of cells to Plot
    tempID  = neurons(eisToPlot(i));  % gets the id number of the ith neuron
    ei = eiFile.getImage(tempID);     % gets ei data for ith neuron, storing the information as a 3D array
    ei = squeeze(ei(1,2:end,:))';
    eiAmps = max(ei)-min(ei);
    [curr_axon_x, curr_axon_y, ~, soma_x, soma_y] = weighted_axon_poly_reg(eiAmps);   
    eiAmps(eiAmps < threshold) = 0;
    
    % banana should this be inside the loop?? LG
    if ~isempty(neuronIDs)
        neuronIDs = unique(neuronIDs);
        data_index = find(unique(neuronIDs) == tempID);
        %data_index = min(data_index);
        
        if ~isempty(data_index)
            soma(data_index, :) = [soma_x soma_y];
            axon_x{data_index, :} = curr_axon_x;
            axon_y{data_index, :} = curr_axon_y;
            max_amp(data_index) = max(eiAmps);
        end
    end
    
    plot(curr_axon_x, curr_axon_y, 'Color', colors(i,:), 'DisplayName', num2str(tempID));
    legend('-DynamicLegend'); hold on;   
    h = scatter(soma_x,soma_y,max(eiAmps),colors(i,:),'filled');
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  
    axis image;
    distances = [];
    for ind = 1:length(curr_axon_x)
        distances = [distances pdist([615 30; curr_axon_x(ind) curr_axon_y(ind)])];
    end
    
    if ~isempty(distances)
        dist_to_elec(i,:) = [tempID min(distances)];
    end 
end
title(sprintf('%0.0f EIs with signal above %0.0f on electrode %0.0f',...
    length(eisToPlot),threshold,electrodeNo)); 

try
    lns = cellfun(@length, axon_x);
catch
    return;
end
x_padded = NaN(numel(axon_x),max(lns));
y_padded = NaN(numel(axon_y),max(lns));

for line = 1:numel(axon_x)
    x_padded(line, 1:lns(line)) = axon_x{line};
    y_padded(line, 1:lns(line)) = axon_y{line};
end

axis off;
disp('Axon IDs and distances to stim electrode');
format shortG
dist_to_elec %#ok<NOPRT>

figure; set(gcf,'Color',0.95*[1 1 1],'Position',[680	510	1050 580]); hold on;
if ~isempty(neuronIDs)
    stim_levels_1 = cell2mat(stim_levels(1,1)); % Get the stim levels for the first neuron
    slmin = 1;
    slmax = size(stim_levels_1,1);  
    vars=struct('IDs', neuronIDs, 'stim', stim_levels, 'activation', ...
        activation,'soma',soma,'axon_x',x_padded,'axon_y', y_padded,...
        'max_amp', max_amp,'electrodeNo',electrodeNo,'positions',positions);    
    hSlider = uicontrol('Style','slider','Min',slmin,'Max',slmax,...
        'SliderStep',[1 1]./(slmax-slmin),'Value',1,...
        'Position',[20 20 200 20], 'UserData', vars); 
    set(hSlider,'Callback',@hsl_callback);
    hListener = addlistener(hSlider,'ContinuousValueChange',@hsl_callback); %#ok<NASGU>
    
    title(strcat('Stimulation amplitude: ',num2str(abs(stim_levels_1(1))), ' uA'));
    scatter(positions(electrodeNo,1),positions(electrodeNo,2),900,'k'); 
    scatter(positions(electrodeNo,1),positions(electrodeNo,2),300,'k', 'filled'); 
    axis image; axis off; 
    plotterfcn(vars,1)
    
    % Add text to identify neuron IDs on the plot
    for n = 1:length(neuronIDs)
        text(double(vars(1).soma(n,1)),double(vars(1).soma(n,2)), num2str(vars(1).IDs(n)));
    end    
  
end

end % End main function

function hsl_callback(hObject,eventdata)
    vars = get(hObject, 'UserData');
    value = get(hObject, 'Value');
    disp(['updating plot,uservalue is ' num2str(value)]);
    plotterfcn(vars, value);    
end
    
function plotterfcn(vars, value)
 cla;
  stimAmpToPlot = vars(1).stim(ceil(value)); 
  colors = jet(20);          % Assign colormap
    for n = 1:length(vars(1).IDs)     % For each cell to be plotted
        stimAmps = vars(n).stim;
        stim_idx = find(stimAmps == stimAmpToPlot,1); % Accounts for different stimAmps for different neurons, adjust to interpolate if needed.     
        if isempty(stim_idx)
            fprintf(['Update code to interpolate for the activation '...
                'level of n%0.0f at stimulation current=%0.4fuA\n'],vars(1).IDs(n),stimAmpToPlot);
        else
            actProb(n) = vars(n).activation(stim_idx); %#ok<AGROW> % Number between 0 and 1           
            colorIdx = ceil(actProb(n)*size(colors,1)); %An index of colors to plot
            if colorIdx == 0; colorIdx = 1; end
            plot(vars(1).axon_x(n,:), vars(1).axon_y(n,:), 'Color', colors(colorIdx,:),'LineWidth',2);
            scatter(vars(1).soma(n,1), vars(1).soma(n,2),vars(1).max_amp(n),colors(colorIdx,:),'filled');
        end
    end
   
    scatter(vars(1).positions(vars(1).electrodeNo,1),vars(1).positions(vars(1).electrodeNo,2),900,'k'); 
    scatter(vars(1).positions(vars(1).electrodeNo,1),vars(1).positions(vars(1).electrodeNo,2),300,'k', 'filled'); 
    c=colorbar('southoutside'); colormap(jet);
    c.Label.String = 'probability of activation at this stimulation level';
    title(['Stimulation amplitude: ',num2str(stimAmpToPlot), ' uA']); 
    axis image; axis off;
    % Banana embed this in the main figure window!
    figure(1000); set(gcf,'Position',[42   111   658   396]); cla; 
    hist(actProb,0.05:0.1:0.95); xlabel('activation probability');ylabel('number of neurons'); 
    title(sprintf('distribution of activation probability at %0.4f uA',stimAmpToPlot)); 
end    
