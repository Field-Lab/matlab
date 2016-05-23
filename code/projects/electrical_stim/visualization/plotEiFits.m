function figure_handle =plotEiFits(pathToEi,varargin)
% PLOTEIFITS(..) plots polynomial fits to the specified EIs. Running
% with only the path to an EI file by default plots all EIs in the
% recording. Can optionally plot a subset of the neurons by passing a list
% of valid neuron IDs to plot.
%
%    Inputs:
%             pathToEi: path to the ei file to use to plot the eis, a string. e.g. '/Volumes/Analysis/2015-10-06-3/data000/data000.ei'
%       Optional:
%         neuronIdList: vector with the neuron IDs to plot. Default is all
%                       neurons in an EI file.
%             somaOnly: true/false, plots only the soma position. default false
%           plotCoords: default true;
%            ei_thresh: default 10;
%       printNeuronIds: prints the IDs of the plotted neurons. default false;
%            plotColor: [];
%            figureNum: [];
%      alternateCoords: can optionally include alternate xy electrode coordinates for plotting over an image, for example default, [];
%
%    Outputs:
%        figure_handle: handle to the figure with plots.
% Lauren Grosberg 5/2016

nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

% Convert ei file into a matrix with the EIs of all the recorded cells
[eiM,neuronIdList_all] = convertEiFileToMatrix(pathToEi);

% Set up defaults for optional parameters
somaOnly = false;
plotCoords = true;
ei_thresh = 10;
printNeuronIds = false;
plotColor = [];
figureNum = [];
alternateCoords = [];
neuronIdList = neuronIdList_all;

% Read the optional input arguments
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'neuronidlist'
            neuronIdList = varargin{j*2};
        case 'somaonly'
            somaOnly = varargin{j*2};
        case 'figurenum'
            figureNum = varargin{j*2};
        case 'plotcoords'
            plotCoords = varargin{j*2};
        case 'ei_thresh'
            ei_thresh = varargin{j*2};
        case 'plotcolor'
            plotColor = varargin{j*2};
        case 'alternatecoords'
            alternateCoords = varargin{j*2};
        case 'printneuronids'
            printNeuronIds = varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

if size(neuronIdList,2) ~=1
    neuronIdList = neuronIdList';
end

% Get electrode coordinates
switch size(eiM,1)
    case 512
        [xc, yc] = getElectrodeCoords512();
        load('adj_mat_512.mat');
        adj_matrix = adj_mat_512;
    case 519
        [xc, yc] = getElectrodeCoords519();
        load('adj_mat_519.mat');
        adj_matrix = adj_mat_519;
        
    case {61,64}
        [xc, yc] = getElectrodeCoords61();
        load('adj_mat_61.mat');
        adj_matrix = adj_mat_61;
    otherwise
        err = MException('MATLAB:InvArgIn',...
            'Unknown array specified - must be 512, 519, or 61');
        throw(err);
end

if ~isempty(alternateCoords)
    dimC = find(size(alternateCoords) == 2);
    if dimC == 1
        xc = alternateCoords(1,:);
        yc = alternateCoords(2,:);
    elseif dimC == 2
        xc = alternateCoords(:,1);
        yc = alternateCoords(:,2);
    end
end

if isempty(figureNum)
    fh = figure;
else
    fh = figure(figureNum);
end

if isempty(plotColor)
    colors = lines(size(neuronIdList,1));
else
    colors = repmat(plotColor,size(neuronIdList,1),1);
end

for n = 1:1:size(neuronIdList,1)
    cellID = neuronIdList(n);
    cellIndex = find(neuronIdList_all == cellID);
    ei = eiM(:,:,cellIndex);  %#ok<FNDSB>
    eiAmps = (max(ei,[],2)-min(ei,[],2))';
    
    [~,col,~] = find(eiAmps > ei_thresh);
    aa = round(eiAmps(col));
    
    sortaa = sort(aa,2,'descend');
    try
        largestAmps = sortaa(1:2);
    catch
        largestAmps = sortaa(1);
    end
    
    [max_amp,max_amp_i] = max(eiAmps);
    adj = adj_matrix{max_amp_i};
    [max_adj,max_adj_i] = max(eiAmps(adj));
    
    indicies = [max_amp_i;adj(max_adj_i)];
    
    % Soma location is center of mass of largest amplitude and largest neighbor
    COMx = sum(xc(indicies).*eiAmps(indicies)) / (max_amp + max_adj);
    COMy = sum(yc(indicies).*eiAmps(indicies)) / (max_amp + max_adj);
    
    valid = 1;
    
    if ~somaOnly
        [XI, YI, ~, COMx, COMy, valid] = weighted_axon_poly_reg(eiAmps);
    end
    
    figure(fh);
    if valid
        if ~somaOnly; plot(XI,YI,'-','Color',colors(n,:)); end
        hold on;
        scatter(COMx,COMy,mean(largestAmps), colors(n,:),'filled');
    else
        fprintf('Did not plot: (above warning for axon %d)\n',cellID);
    end
    if printNeuronIds
        text(double(COMx),double(COMy),num2str(cellID));
    end
    %     hold on; scatter(yc(row),xc(row),eiAmps(row)*6,colors(n,:),'filled');   % Plot eis
    
end
axis image; axis off;
if plotCoords
    hold on; scatter(xc,yc,5,'black','filled');
end
figure_handle = figure(fh);
end % End function.
