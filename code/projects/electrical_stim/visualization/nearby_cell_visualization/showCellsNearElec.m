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
    
    % Plotting ..
    % Overlay rf fits and IDs
    figure; set(gcf,'Color',[1 1 1]); hold on;
    
 
    plot_rf_summaries(datarun, ONPcellIDs, 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b');
    plot_rf_summaries(datarun, OFFPcellIDs, 'clear', false, 'label', true, 'label_color', 'k', 'plot_fits', true, 'fit_color', 'k');
    plot_rf_summaries(datarun, ONMcellIDs, 'clear', false, 'label', true, 'label_color', 'm', 'plot_fits', true, 'fit_color', 'm');
    plot_rf_summaries(datarun, OFFMcellIDs, 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'g');
    axis image; axis off;
    % legend('ON parasol','OFF parasol','ON midget','OFF midget')
    title(sprintf('ON parasol (blue) OFF parasol (black)\n ON midget (pink) OFF midget (green)\n stimulating electrode %0.0f',electrodeNo));
    %{
    figure
    for n = 1:4
        switch n
            case 1 
                plot_cells = ONPcellIDs;
                label_color = 'b';
            case 2
                plot_cells = OFFPcellIDs;
                label_color = 'k';
            case 3
                plot_cells = ONMcellIDs;
                label_color = 'm';
            case 4
                plot_cells = OFFMcellIDs;
                label_color = 'g';
            otherwise    
        end
        for c = 1:length(plot_cells)
            ei = eiFile.getImage(plot_cells(c));  % gets ei data for ith neuron, storing the information as a 3D array
            ei = squeeze(ei(1,2:end,:))';
            eiAmps = max(ei)-min(ei); 
            [axon_x, axon_y, soma_x, soma_y] = weighted_axon_poly_reg(eiAmps);
            plot(axon_x, axon_y, label_color); hold on;
            scatter(soma_x,soma_y,max(eiAmps), label_color,'filled'); hold on;
            text(double(soma_x),double(soma_y),num2str(plot_cells(c))); hold on;
        end
    end  
    %} 
catch [cellSpec, cellType] = get_cell_indices( datarun, 'all' );
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
    % Plotting ..
    % Overlay rf fits and IDs
    figure; set(gcf,'Color',[1 1 1]); hold on;
    plot_rf_summaries(datarun, cellIDs, 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b');
    axis image; axis off;
    % legend('ON parasol','OFF parasol','ON midget','OFF midget')
    title(sprintf('%s stimulating electrode %0.0f',datarun.names.short_name,electrodeNo));
    
end %end catch statement
figure; set(gcf,'Color',[1 1 1]); hold on;
eisToPlot = get_cell_indices(datarun,cellIDs);
neurons = datarun.cell_ids(:)';
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
    axis image;
end

text(positions(maxRecElec(cellIDs),1),positions(maxRecElec(cellIDs),2),num2str(cellIDs),'HorizontalAlignment','center');
scatter(positions(electrodeNo,1),positions(electrodeNo,2),50,'k'); axis off;
text(positions(electrodeNo,1),positions(electrodeNo,2),'stim elec');


figure; set(gcf,'Color',[1 1 1]); hold on;

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
            load(strcat(dataPath, fileName));
            if autoFile
                neuronID = elecRespAuto.neuronInfo.neuronIds;
            else    
                neuronID = elecResp.cells.main;
            end
            if ismember(neuronID, cellIDs)
                neuronIDs(k) = neuronID;
                if autoFile
                    stim_levels(k, :) = elecRespAuto.stimInfo.listAmps;
                    activation(k, :) = elecRespAuto.LogisticReg;
                else    
                    stim_levels(k, :) = elecResp.stimInfo.stimAmps;
                    activation(k, :) = elecResp.analysis.successRates;
                end
                k = k + 1;
            end
        end
    end
end

for i = 1:length(eisToPlot)
    tempID  = neurons(eisToPlot(i));          % gets the id number of the ith neuron
    ei = eiFile.getImage(tempID);  % gets ei data for ith neuron, storing the information as a 3D array
    ei = squeeze(ei(1,2:end,:))';
    eiAmps = max(ei)-min(ei);
    [curr_axon_x, curr_axon_y, ~, soma_x, soma_y] = weighted_axon_poly_reg(eiAmps);
    
    eiAmps(find(eiAmps<threshold)) = 0;
    
    %{
    if isempty(curr_axon_x)
        curr_axon_x = soma_x; %Something to plot, so it shows up in the legend
        curr_axon_y = soma_y;
    end
    %}
    
    if ~isempty(neuronIDs)
        %neuronIDs = unique(neuronIDs);
        data_index = find(neuronIDs == tempID);
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
dist_to_elec

figure; set(gcf,'Color',[1 1 1]); hold on;
if ~isempty(neuronIDs)

    slmin = 1;
    slmax = length(stim_levels);
    
    vars=struct('IDs', neuronIDs, 'stim', stim_levels, 'activation', activation, ...
        'soma', soma, 'axon_x', x_padded, 'axon_y', y_padded, 'max_amp', max_amp);
    
    hSlider = uicontrol('Style','slider','Min',slmin,'Max',slmax,...
        'SliderStep',[1 1]./(slmax-slmin),'Value',1,...
        'Position',[20 20 200 20], 'UserData', vars);
    

    
    set(hSlider,'Callback',@hsl_callback);
    hListener = addlistener(hSlider,'ContinuousValueChange',@hsl_callback);
    title(strcat('Stimulation amplitude: ',num2str(abs(stim_levels(1,1))), ' uA'));
    scatter(positions(electrodeNo,1),positions(electrodeNo,2),900,'k'); 
    scatter(positions(electrodeNo,1),positions(electrodeNo,2),300,'k', 'filled'); 
    plotterfcn(vars,1)
    
    for n = 1:length(neuronIDs)
        text(double(vars.soma(n,1)),double(vars.soma(n,2)), num2str(vars.IDs(n)));
    end    

   
end

end

function hsl_callback(hObject,eventdata)
    vars = get(hObject, 'UserData');
    value = get(hObject, 'Value');
    plotterfcn(vars, value);
end
    
function plotterfcn(vars, value)
  stim_int = size(vars.stim,2);
  colors = gray(stim_int);
    for n = 1:size(vars.stim, 1)
        intensity = vars.activation(n,round(value));
        gray_ind = stim_int - (ceil(intensity*(stim_int - 1)));
        plot(vars.axon_x(n,:), vars.axon_y(n,:), 'Color', colors(gray_ind,:));
   
        scatter(vars.soma(n,1), vars.soma(n,2),vars.max_amp(n),colors(gray_ind,:),'filled');   
    end
    title(strcat('Stimulation amplitude: ',num2str(abs(vars.stim(n,round(value)))), ' uA'));
end    
