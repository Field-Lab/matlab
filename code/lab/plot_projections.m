function image = plot_projections(projections, varargin)
%image = PLOT_PROJECTIONS(projections, [params]) 
%   plots projections in the current figure or axes. If clusters are
%   specified, they will be plotted on top of the projections.
%
% usage:  axes = plot_projections(projections, 'clusters', 'data000.model', 'electrode', 1)
%
% arguments:   projections  -- a string specifiying the path to a
%              projections file, or a matrix of size N_DIMENSIONS x
%              N_POINTS.
%
% optional arguments:
%              clusters     -- a string specifying a path to a model file,
%              or an EMModel Java object from a single electrode obtained
%              by calling: 
%
%                 model = edu.ucsc.neurobiology.vision.io.ClusteringModelFile('path.model')
%                 model.getNeuronExtraction(electrode)
%
%              electrode -- if paths are provided to projections and/or
%              model files, an electrode to view must be specified
%
%              points -- number of points to plot in the PC projections
%              under the clusters. the default is 1000 points.
%
%
% tamachado@salk.edu 6/22/09

% parse parameters
p = inputParser;
EMPTY = -1;

p.addParamValue('points', 1000); 
p.addParamValue('clusters', EMPTY); 
p.addParamValue('electrode', EMPTY); 

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% check if electrode has been specified
if isa(params.clusters,'char') || isa(projections,'char')
    if params.electrode == EMPTY
       error('If you specify file paths, you must specify an electrode.')
    end
end

if params.clusters ~= EMPTY
    % load the clusters if necessary
    if isa(params.clusters,'char')
        model = edu.ucsc.neurobiology.vision.io.ClusteringModelFile(params.clusters);
        params.clusters = model.getNeuronExtraction(params.electrode);
    end
    cRow = 2;
else
    % make one row of subplots if there are no clusters to plot
    cRow = 1;
end

if isa(projections,'char')
    % open the projections file
    prjObject = edu.ucsc.neurobiology.vision.matlab.ReadProjections(projections);
    % get the projections from the specified electrode
    prjObject.readElectrode(params.electrode);
    projections = prjObject.getProjections();
    % get the number of spikes on that electrode
    nSpikes = prjObject.getSpikeCount();
    % remove trailing zeros
    projections = projections(:,1:nSpikes);
end

% Plot projections
density_plot_axes = subplot_axes(gcf,[0 0 1 1],0.1,0.2,3,cRow);

% Generate the density plot
spike_proj = projections';
dim_to_plot = [1 2;1 3;2 3];
bins = 230;
dimension_prefix = 'PC';
if params.points > length(spike_proj), params.points = length(spike_proj); end

density_gradient = cell(1, length(dim_to_plot));
hist_centers = cell(1, length(dim_to_plot));
density_gradient_norm_factor = cell(1, length(dim_to_plot));

for ff=1:length(dim_to_plot)
    % note which dimensions to plot (e.g. PC 1 and PC 2)
    first_dimension = dim_to_plot(ff,1);
    second_dimension = dim_to_plot(ff,2);
    dim = [first_dimension second_dimension];
    
    %generate density gradient, and save histogram centers
    %make 2d histogram of spike_locations
    [temp,hist_centers{ff}]=hist3([spike_proj(:,first_dimension) spike_proj(:,second_dimension)],[1 1]*bins);

    % save histogram maximum value
    density_gradient_norm_factor{ff} = max(max(temp));
    
    % find where that location is
    [x,y] = find(temp == density_gradient_norm_factor{ff});
  
    %normalize the values
    temp=temp/density_gradient_norm_factor{ff};
    density_gradient{ff}=-1*((log(0.2*temp+0.01)+2*temp+6).*(temp>0));
    density_gradient{ff}=rot90(density_gradient{ff});

    %plot the density gradient
    d = density_gradient{ff};
    axes(density_plot_axes{ff});
    h = imagesc(d); axis square;
    
    % label axes
    xlabel(sprintf('%s %d',dimension_prefix,first_dimension));
    ylabel(sprintf('%s %d',dimension_prefix,second_dimension));
    colormap('bone');
    
    %remove tick marks
    set(density_plot_axes{ff},'YTick',[],'XTick',[]);

    % plot clusters if desired
    if params.clusters ~= EMPTY
        axes(density_plot_axes{ff+3});
        plot(spike_proj(1:params.points,dim(1)),spike_proj(1:params.points,dim(2)),'m.'); hold on;
        rng = [xlim ylim];
        nClusters = length(params.clusters.probability);
        for cc = 1:nClusters
            % get the pdf for the cluster
            cm = params.clusters.means(cc,dim);
            gg = gmdistribution(cm,params.clusters.covariances(cc,dim));
            fn = @(x,y) pdf(gg,[x y])./max(pdf(gg,[x y]));
            % plot it
            ezcontour(fn,rng,params.points);
        end
        title(sprintf('nClusters = %d',nClusters))
        axis square;
        %remove tick marks
        set(gca,'YTick',[],'XTick',[]);
        
        % label axes
        xlabel(sprintf('%s %d',dimension_prefix,first_dimension));
        ylabel(sprintf('%s %d',dimension_prefix,second_dimension));
    end

end

% save values in a structure
image.density_gradient = density_gradient;
image.density_gradient_norm_factor = density_gradient_norm_factor;
image.histogram_centers = hist_centers;
