function fittedGLM=glm_fit_from_WNrun(cells, dataset, stim_description, stim_length, d_save)
% glm_fit_from_classrun(cells, dataset, stim_description, optional: stim_length, d_save)
% dataset='2014-11-05-2/data009_nps';
% cells={2372,2523}
% stim_description 'RGB-10-2-0.48-11111-32x32'
% stim_length, optional, default 30 min
% dsave, optional, where to save the fit

if nargin==3
    stim_length=1800; 
    d_save='/Volumes/Analysis/nora/colorglmfits';
elseif nargin==4
    d_save='/Volumes/Analysis/nora/colorglmfits';
end

%% Stimulus and GLM parameters

% Get stimulus parameters from descriptor
xml_file=['/Volumes/Analysis/stimuli/white-noise-xml/' stim_description '.xml'];
dashes=find(stim_description=='-');
StimulusPars.type=stim_description(1:dashes(1)-1);
StimulusPars.pixelsize = str2double(stim_description(dashes(1)+1:dashes(2)-1));
StimulusPars.refreshrate = str2double(stim_description(dashes(2)+1:dashes(3)-1));
StimulusPars.RNG = str2double(stim_description(dashes(4)+1:dashes(5)-1));
try
    StimulusPars.height = str2double(stim_description(dashes(5)+1:dashes(5)+2)); 
    StimulusPars.width = str2double(stim_description(dashes(5)+4:dashes(5)+5));
catch
    StimulusPars.height = 32; 
    StimulusPars.width = 64;
end
StimulusPars.tstim = 1/120;
fitframes=stim_length*120/StimulusPars.refreshrate; % seconds * 120 frames per second / interval

% Load a bunch of parameters about the GLM fit type, default is rk1 GLM fit
GLM_pars

% Load datarun
datarun=load_data(dataset);
datarun=load_neurons(datarun);
datarun=load_sta(datarun);
datarun=load_params(datarun);

% Make the directory to save
if ~exist(d_save,'dir'), mkdir(d_save); end
GLMType.d_save = d_save;

%% Load Movie
disp('Loading Stimulus Movies')
[temp_fitmovie,height,width,~,~] = get_movie(xml_file, datarun.triggers, fitframes/2);
temp_fitmovie=permute(temp_fitmovie,[2 1 3 4]);
fitmovie_color=zeros(width,height,3,fitframes);
for i=1:fitframes
    fitmovie_color(:,:,:,i)=temp_fitmovie(:,:,:,ceil(i/2));
end
clear temp_fitmovie height width i
% testmovie = get_rawmovie(raw_file, testframes);

%% Load Cell Specific Elements   Spikes and STA
for i_cell = 1:length(cells)
    clear glm_cellstruct
    cid = cells{i_cell};
    %[celltype , cell_savename, ~]  = findcelltype(cid, datarun.cell_types);
    
    if ~exist(sprintf('%s/%s.mat', d_save,num2str(cid)),'file')
        
        % Load cell info
        glm_cellinfo.cid           = cid;
        %glm_cellinfo.exp_nm        = exp_nm;
        %glm_cellinfo.celltype      = celltype;
        glm_cellinfo.cell_savename = num2str(cid);
        glm_cellinfo.fitname       = GLMType.fitname;
        glm_cellinfo.computedtstim = StimulusPars.tstim;
        master_idx         = find(datarun.cell_ids == cid);
        stafit_centercoord = ( datarun.vision.sta_fits{master_idx}.mean );
        stafit_sd          = ( datarun.vision.sta_fits{master_idx}.sd   );
        slvdim.height      = StimulusPars.height; slvdim.width = StimulusPars.width;
        [center_coord,~]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, StimulusPars, slvdim);
        glm_cellinfo.slave_centercoord = center_coord;
        RGB=RGB_weights(datarun,master_idx);
        glm_cellinfo.RGB=RGB;
        clear stafit_centercoord slvdim sd
        
        % Turn RGB movie into greyscale movie
        fitmovie=squeeze(RGB(1)*fitmovie_color(:,:,1,:)+ ...
            RGB(2)*fitmovie_color(:,:,2,:)+ ...
            RGB(3)*fitmovie_color(:,:,3,:));
        clear RGB
        
        % Spike loading
        spikes=datarun.spikes{master_idx};
        glm_cellinfo.WN_STA = datarun.stas.stas{master_idx};
        clear cell_savename
        
        % Align the spikes and the movies;
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
        spike.home=spikes_adj;
        clear spikes_adj;
        
        % Execute GLM
        tic
        [fittedGLM]     = glm_execute_CP(GLMType, spike,0, fitmovie, glm_cellinfo);
        toc
        
        %{
        xvalperformance = eval_xvalperformance_NEW_CP(fittedGLM, StimulusPars.slv, cell_organizedspikes,neighbor_organizedspikes,testmovie);
        fittedGLM.xvalperformance  = xvalperformance;
        %}
        fittedGLM.d_save           = d_save;
        eval(sprintf('save %s/%s.mat fittedGLM', d_save, glm_cellinfo.cell_savename));
        
    else
        error('Previous results still in directory')
        
    end
    
end
