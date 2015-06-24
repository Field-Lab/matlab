%% NB 2015-05-01 Fitting the GLM from a white noise run

%
% PATHS NEEDED
% Vision.jar, such as javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar')
% The lab codebase, addpath(genpath('Repo location /matlab/code/lab'))
% The glm code folder, addpath(genpath('Repo location /matlab/code/projects/glm))
%
% TROUBLESHOOTING
% If it is not working, try using 'center_verification', true to set the
%   cell's location yourself.

function [fittedGLM] = glm_fit_from_WN(cells, dataset, stim, varargin)
% DESCRIPTION
% This is a wrap that loads the white noise movie specified by the XML file
% given by the stim_description. It also loads the spikes and STA
% information from the datarun. Then it passes these things to the glm_fit functions.
% If you already have the movie loaded up, you can specify stim as that instead of a
% string. If you give an RGB movie, it will be converted to BW by taking RGB
% weights from the STA and making a psuedoBW movie.
%
% The GLM architecture and settings can be changed in
% glm_parameters.m
%
% USAGE
% glm_fit_from_WN(cells, dataset, stim, optional: stim_length, d_save)
% example
% fittedGLM = glm_fit_from_WN([121], '2012-08-09-3/data002','RGB-8-1-0.48-11111')

% INPUTS

%  REQUIRED

%   cells
%       a vector of the cell ids to fit. If you do this, you should set
%       dsave so that the fits save to disk
%       ex cells=[2372,2523]
%   dataset
%       a string with dataset name. exactly what you would put into the
%       load_data function
%       ex ataset ='2014-11-05-2/data009_nps';
%   stim: either a string like 'RGB-10-2-0.48-11111-32x32'
%       or the matfile movie
%       with a frame for every 1/120 seconds 
%       (regardless of interval, so will have repeated frames if interval >1)
%       and stim size x time if black and white or stim size x 3 x time if RGB

% OPTIONAL KEYWORDS
% stim_length, optional, seconds, default is 900
% dsave, optional, where to save the fit. If no path is given, the fit will
%   not be saved to disk
% monitor_refresh: usually should be 120Hz. This is NOT the interval!
%   Just the monitor speed!
% center_verification: default false
%   set to TRUE to pick the cell location youself, instead of relying on
%   the automatic center picking. A figure will come up with instructions
%   of when and where to click.

% Parse optional input
p = inputParser;
p.addParameter('stim_length', 900)
p.addParameter('d_save', 0)
p.addParameter('monitor_refresh', 120)
p.addParameter('center_verification', true)
p.parse(varargin{:});
stim_length = p.Results.stim_length;
d_save = p.Results.d_save;
monitor_refresh = p.Results.monitor_refresh;
center_verification = p.Results.center_verification;
clear p
frames_per_trigger = 100;

% Load datarun
datarun=load_data(dataset);
datarun=load_neurons(datarun);
datarun=load_sta(datarun, 'load_sta', cells);
datarun=load_params(datarun);

% Make the directory to save
if isstr(d_save)
    if ~exist(d_save,'dir'), mkdir(d_save); end
end

% Get Movie
if isstr(stim)
    disp(['Using ' stim ' XML file'])
    xml_file=['/Volumes/Analysis/stimuli/white-noise-xml/' stim '.xml'];
    dashes=find(stim=='-');
    WN_type=stim(1:dashes(1)-1);
    interval = str2double(stim(dashes(2)+1:dashes(3)-1));
    fitframes=stim_length*monitor_refresh/interval;
    disp('Loading Stimulus Movies')
    [temp_fitmovie,height,width,~,~] = get_movie(xml_file, datarun.triggers, fitframes);
    fitmovie_color=zeros(height,width,3,fitframes*interval);
    for i=1:(fitframes*interval)
        fitmovie_color(:,:,:,i)=temp_fitmovie(:,:,:,ceil(i/interval));
    end
    if strcmp(WN_type, 'BW')
        fitmovie = squeeze(fitmovie_color(:,:,1,:));
    end
else
    fitmovie = stim;
end
clear temp_fitmovie height width i

%% Load Cell Specific Elements Spikes and STA

for i_cell = 1:length(cells)
    cid = cells(i_cell);
    
    % Cell Info
    glm_cellinfo.cid           = cid;
    glm_cellinfo.cell_savename = num2str(cid);
    master_idx         = find(datarun.cell_ids == cid);
   
    % Make the movie psuedo BW if it isn't already BW
    if ~exist('fitmovie', 'var')
        [RGB, fitmovie] = RGB_to_BW(datarun, master_idx, 'color_movie', fitmovie_color);
        clear fitmovie_color
        fitmovie = permute(fitmovie, [2 1 3]);
    end
    
    % Check the size of the movie
    size_fitmovie = size(fitmovie);
    assert(length(size_fitmovie) == 3, ['The movie size is ' num2str(size_fitmovie)])
    
    % Spike loading
    spikes = datarun.spikes{master_idx};
    
    % Load the STA and check the size
    WN_STA = datarun.stas.stas{master_idx};
    if size(WN_STA,3) == 3
        temp_STA=squeeze(RGB(1)*WN_STA(:,:,1,:)+ ...
            RGB(2)*WN_STA(:,:,2,:)+ ...
            RGB(3)*WN_STA(:,:,3,:));
        clear WN_STA
        WN_STA = permute(temp_STA, [2 1 3]);
        clear temp_STA
    end
    WN_STA = squeeze(WN_STA);
    size_STA = size(WN_STA);
    assert(length(size_STA) == 3, ['The STA size is ' num2str(size_STA)])  
    if size_fitmovie(1:2) ~= size_STA(1:2),
        warn('Fitmovie and STA are not the same size')
    end
    
    % Pull out the cell location
    center_coord = datarun.vision.sta_fits{master_idx}.mean;
    center(1) = round(center_coord(1)); % x_coord
    center(2) = size(fitmovie,2) - round(center_coord(2)); %y_coord
    clear cell_savename
    
    % Align the spikes and the movies to the triggers;
    spikes_adj=spikes;
    n_block=0;
    for i=1:(length(datarun.triggers)-1)
        actual_t_start=datarun.triggers(i);
        supposed_t_start=n_block*frames_per_trigger/monitor_refresh;
        idx1=spikes > actual_t_start;
        idx2=spikes < datarun.triggers(i+1);
        spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
        n_block=n_block+1;
    end
    clear spikes
    fitspikes=spikes_adj;
    clear spikes_adj;
    
    % If coupling, load up the neighboring spikes and adjust them to the
    % triggers
    [GLMT, ~] = glm_parameters;
    if GLMT.CouplingFilters
        cell_ids = zeros(2, length(datarun.spikes));
        cell_ids(2, :) = 1:length(datarun.spikes);
        cell_ids(1, :) = datarun.cell_ids;
        % if you want to only pair a certain cell type, you'd have to change
        % this to only include the cell_ids and sta_fits of the candidate
        % cells 
        paired_cells = subR_pick_neighbor_cells(center_coord, cell_ids, datarun.vision.sta_fits);
        for i_pair = 1:length(paired_cells)
            n_block=0;
            spikes = datarun.spikes{datarun.cell_ids == paired_cells(i_pair)}; 
            for i=1:(length(datarun.triggers)-1)
                actual_t_start=datarun.triggers(i);
                supposed_t_start=n_block*frames_per_trigger/monitor_refresh;
                idx1=spikes > actual_t_start;
                idx2=spikes < datarun.triggers(i+1);
                spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
                n_block=n_block+1;
            end
            neighborspikes{i_pair} = spikes_adj;
        end
    else
        neighborspikes = 0;
    end
    
    % Check that everything is working, and automatically find the center
    [STA, center] = STA_Test(fitspikes, fitmovie, center_verification);
    
    % Execute and save GLM
    tic
    fittedGLM     = glm_fit(fitspikes, fitmovie, center, 'WN_STA', WN_STA, 'monitor_refresh', monitor_refresh, 'neighborspikes', neighborspikes);
    toc
    
    if GLMT.CouplingFilters
        fittedGLM.paired_cells = paired_cells;
    end
    if isstr(d_save)
        eval(sprintf('save %s/%s.mat fittedGLM', d_save, glm_cellinfo.cell_savename));
    end
    
    
end


end

function paired_cells=subR_pick_neighbor_cells(mean, cell_ids, sta_fits)
    
     [~, GLMPars] = glm_parameters;
     NumCells = length(cell_ids);
     distance=zeros(NumCells,1);
     
     % Calculate distance between RFs
     for i_pair=1:NumCells
         distance(i_pair)=norm(sta_fits{cell_ids(2,i_pair),1}.mean-mean);
         if distance(i_pair)==0
             distance(i_pair)=NaN; % if it is the cell itself, set distance to infinity
         end
     end
     
     % Choose the closest cells
     [~,indices]=sort(distance);
     paired_cells=cell_ids(1,indices(1:GLMPars.spikefilters.cp.n_couplings));

end
