clear; close all;  clc

dataset='2014-11-24-3/data005';

cells={96};

d_save='/Volumes/Analysis/nora/colorglmfits_test';
xml_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
fitframes=20*60*120/2; % 30 minutes * 60 seconds * 120 frames per second / interval of 2

%% DICTATE GLMTYPE and Datasets and cells  EDITS DONE HERE! 
GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
%GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
%GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; GLMType.fixedSPlength = 13;  GLMType.fixedSP_nullpoint = 'mean'; 
GLMType.nullpoint = 'mean';
GLMType.fit_type = 'WN'; GLMType.map_type = 'mapPRJ';
GLMType.debug = false;
GLMType.specialchange = false;
GLMType.CBP=false;
GLMType.stimfilter_mode = 'rk1';
%GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.input_pt_nonlinearity      = false;
%GLMType.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
GLMType.input_pt_nonlinearity_type = 'raisepower_meanafter';
GLMType.CONVEX = false;
GLMType.DoubleOpt = false;
GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.Subunits = false;
GLMType.Saccades=false;
GLMType.func_sname = 'glmwrap24_CP';
GLMType.fullmfilename =mfilename('fullpath');
i_exp = 1; i_cell = 1;
GLMType.fitname  = GLM_fitname(GLMType);

% Load dataruns
datarun=load_data(dataset);
datarun=load_neurons(datarun);
datarun=load_sta(datarun);
datarun=load_params(datarun);

% Stimulus details
StimulusPars.pixelsize = 10;
StimulusPars.height = 32; StimulusPars.width  = 64;
StimulusPars.refreshrate = 2;
%StimulusPars.frames_pertrigger = 50;
StimulusPars.tstim = 1/120;
StimulusPars.type = 'RGB';
StimulusPars.RNG  = 11111;

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
        
        %% STA test
        STA=zeros(size(fitmovie,1),size(fitmovie,2),30);
        for i=1:length(spike.home)
            spike_fr=ceil(spike.home(i)*120);
            if spike_fr > 29 && spike_fr < size(fitmovie,3)
            STA=STA+fitmovie(:,:,(spike_fr-29):spike_fr);
            end
        end
        
        %%
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
