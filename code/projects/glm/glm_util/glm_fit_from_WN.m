%% NB 2015-05-01 Fitting the GLM from a white noise run
% This is a wrap that loads the white noise movie specified by the XML file
% given by the stim_description. It also loads the spikes and STA
% information from the datarun. Then it passes these things to glm_fit. 
%
% The GLM architecture and settings can be changed in
% glm_parameters.m
%
% example fittedGLM = glm_fit_from_WN([121], '2012-08-09-3/data002',
% 'RGB-8-1-0.48-11111')
%
% PATHS NEEDED
% Vision.jar such as javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar')
% The lab codebase, addpath(genpath('Repo location /matlab/code/lab'))
% The glm code folder, addpath(genpath('Repo location /matlab/code/projects/glm))

function [fittedGLM] = glm_fit_from_WN(cells, dataset, stim_description, varargin)

% INPUTS

% REQUIRED
% glm_fit_from_WN(cells, dataset, stim_description, optional: stim_length, d_save)
% dataset='2014-11-05-2/data009_nps';
% cells=[2372,2523]
% stim_description 'RGB-10-2-0.48-11111-32x32'

% OPTIONAL KEYWORDS
% stim_length, optional, default 15 min
% dsave, optional, where to save the fit. If no path is given, the fit will
% not be saved to disk
%

% Parse optional input
p = inputParser;
p.addParamValue('stim_length', 900)
p.addParamValue('d_save', 0)
p.parse(varargin{:});
stim_length = p.Results.stim_length;
d_save = p.Results.d_save;
clear p

% Load datarun
datarun=load_data(dataset);
datarun=load_neurons(datarun);
datarun=load_sta(datarun, 'load_sta', cells);
datarun=load_params(datarun);

% Make the directory to save
if isstr(d_save)
    if ~exist(d_save,'dir'), mkdir(d_save); end
end

% Get Movie Params
disp(['Using ' stim_description ' XML file'])
xml_file=['/Volumes/Analysis/stimuli/white-noise-xml/' stim_description '.xml'];
dashes=find(stim_description=='-');
WN_type=stim_description(1:dashes(1)-1);
interval = str2double(stim_description(dashes(2)+1:dashes(3)-1));
fitframes=stim_length*120/interval;

% Load Movie
disp('Loading Stimulus Movies')
[temp_fitmovie,height,width,~,~] = get_movie(xml_file, datarun.triggers, fitframes);
fitmovie_color=zeros(height,width,3,fitframes);
for i=1:fitframes
    fitmovie_color(:,:,:,i)=temp_fitmovie(:,:,:,ceil(i/interval));
end
disp(size(fitmovie_color))
clear temp_fitmovie height width i

%% Load Cell Specific Elements   Spikes and STA
for i_cell = 1:length(cells)
    cid = cells(i_cell);
    
    % Cell Info
    glm_cellinfo.cid           = cid;
    glm_cellinfo.cell_savename = num2str(cid);
    master_idx         = find(datarun.cell_ids == cid);
    
    % RGB weights
    RGB = RGB_weights(datarun,master_idx);
    glm_cellinfo.RGB=RGB;
    
    % Turn RGB movie into greyscale movie
    fitmovie=squeeze(RGB(1)*fitmovie_color(:,:,1,:)+ ...
        RGB(2)*fitmovie_color(:,:,2,:)+ ...
        RGB(3)*fitmovie_color(:,:,3,:));
    
    % Spike loading
    spikes = datarun.spikes{master_idx};
    WN_STA = datarun.stas.stas{master_idx};
    if size(WN_STA,3) == 3
        temp_STA=squeeze(RGB(1)*WN_STA(:,:,1,:)+ ...
            RGB(2)*WN_STA(:,:,2,:)+ ...
            RGB(3)*WN_STA(:,:,3,:));
        clear WN_STA
        WN_STA = temp_STA;
        clear temp_STA
    end
    center_coord = round(datarun.vision.sta_fits{master_idx}.mean);
    center.y_coord = center_coord(1);
    center.x_coord = size(fitmovie,1) - center_coord(2);
    clear cell_savename
    
    % Align the spikes and the movies to the triggers;
    spikes_adj=spikes;
    n_block=0;
    for i=1:(length(datarun.triggers)-1)
        actual_t_start=datarun.triggers(i);
        supposed_t_start=n_block*100/120;
        idx1=spikes > actual_t_start;
        idx2=spikes < datarun.triggers(i+1);
        spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
        n_block=n_block+1;
    end
    clear spikes
    fitspikes=spikes_adj;
    clear spikes_adj;
    
    % Execute and save GLM
    tic
        fittedGLM     = glm_fit(fitspikes, fitmovie, center, 'WN_STA', WN_STA);
    toc
    if isstr(d_save)
        eval(sprintf('save %s/%s.mat fittedGLM', d_save, glm_cellinfo.cell_savename));
    end
    
    
end
