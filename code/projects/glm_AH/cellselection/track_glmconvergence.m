% AKHEITMAN 2014-11-13
% Rewriting original code to get down the glm_convergence
% Loads up GLM BPS scores for various percentages of fit data
% Keeps both the raw log prob as well as the bps at each convergence point


clear; close all;  clc
% Use this for the convergence check!
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS


% SETUP cells and experiments, the TYPE of GLM (GLMType) 
%debugit = true;
BD = NSEM_BaseDirectories;
exptests = [1 2 3 4];
pct_usage = [.05:.1:.95];
pct_usage = [pct_usage 1];

glmconvdir = sprintf('%s/glm_convergence', BD.Cell_Selection);
if ~exist(glmconvdir,'dir'), mkdir(glmconvdir); end 
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));

%%
glm_conv = allcells;

GLMType.debug = false;
GLMType.fit_type = 'NSEM'; GLMType.map_type = 'mapPRJ'; 
GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
GLMType.nullpoint = 'mean'; 
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.CONVEX = true; % with relation to the filters .. are parameters used linearly in the GLM. 
GLMType.specialchange = true;
GLMType.specialchange_name = 'Fit_Convergence';
GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
GLMType.func_sname = 'glmwrap_24';
GLMType.fullmfilename =mfilename('fullpath'); 
GLMType.fitname0  = GLM_fitname(GLMType);
i_exp = 1; i_cell = 1; i_pct = .05;
convergence_test = cell(4,1);
for i_exp = 1:4
    glm_conv{i_exp}.scores_stimtype_ind = 'One is White Noise, Two is NSEM';
    glm_conv{i_exp}.scores_celltype_ind = 'One is On Parasols, Two is Off Parasols';
    glm_conv{i_exp}.normedbps_note  = 'Rows are Cells, Columns are Percent of fit data used, Normed BPS';
    glm_conv{i_exp}.fit_pct         = pct_usage;
    glm_conv{i_exp}.fitname         = GLMType.fitname0;
end

i_exp = 1; i_stimtype = 1; i_celltype = 1; 

%% Setup looping structure % Collect Raster Scores
for i_exp = 1:4    
    for i_stimtype = 1:2    
        for i_celltype = 1:2
            exp_nm  = allcells{i_exp}.exp_nm;
            expname = allcells{i_exp}.expname;
            if i_stimtype == 1, stimtype = 'WN';   end
            if i_stimtype == 2, stimtype = 'NSEM'; end
            inputs.exp_nm       = exp_nm; 
            inputs.map_type     = 'mapPRJ';
            inputs.stim_type    = stimtype;
            
            if i_celltype == 1; cells = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cells = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
            
            if exist('debugit','var') && debugit
                cells = cells(1:3); 
               % pct_usage = pct_usage((end-3):end);
            end
            glm_conv{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.normedbps   = zeros(length(cells),length(pct_usage));
            glm_conv{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.raw_logprob = zeros(length(cells),length(pct_usage));
            
            for i_pct = 1:length(pct_usage)
                pct = pct_usage(i_pct);
                percentile = round(100*pct);
                inputs.fitname      = sprintf('%s_%dPct',GLMType.fitname0,percentile);
                load_dir = NSEM_secondaryDirectories('savedir_GLMfit', inputs);
                display(sprintf('Working on %s %s Percent Fit %d', stimtype,exp_nm , percentile ) )
                
                for i_cell = 1:length(cells)
                    cid = cells(i_cell);
                    cell_savename = sprintf('%s_%d', celltype,cid);
                    %display(sprintf('Working on %s  Cell: %s', exp_nm,cell_savename));
                    eval(sprintf('load %s/%s.mat fittedGLM', load_dir, cell_savename));
                    convergence_test{i_exp}.normedbps(i_cell,i_pct) = fittedGLM.xvalperformance.glm_normedbits;
                    
                    glm_conv{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.normedbps(i_cell,i_pct)       = fittedGLM.xvalperformance.glm_normedbits;
                    glm_conv{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.raw_logprob(i_cell,i_pct) = fittedGLM.xvalperformance.logprob_glm_raw;
                end
            end
        end
    end
end

eval(sprintf('save %s/glm_conv.mat glm_conv', glmconvdir))

%% Normalize the Matrix  so that we have usable percentages
glm_conv2 = glm_conv;
for i_exp = 1:4
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    
    cell_savename = sprintf('ONPar_%d', allcells{i_exp}.ONP(1)); 
    for i_stimtype = 1:2
        if i_stimtype == 1, stimtype = 'WN';   end
        if i_stimtype == 2, stimtype = 'NSEM'; end
        [StimPars] = Directories_Params_v23(exp_nm, stimtype, 'mapPRJ');
        
        glm_conv2{i_exp}.scores.stim_type{i_stimtype}.framecount = length(StimPars.slv.testframes);
    end
end

for i_exp = 1:4
    for i_stimtype = 1:2
        framecount = glm_conv2{i_exp}.scores.stim_type{i_stimtype}.framecount
        for i_celltype = 1:2
            
            scores = glm_conv2{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype};
            for i_metric = 1:2
                if i_metric == 1, a_mat = scores.normedbps;     end
                if i_metric == 2, a_mat = exp((scores.raw_logprob/framecount));  scores.probperframe = a_mat; end
                
                
                
                cellnum = size(a_mat,1);
                for i_row = 1:cellnum          
                	a_mat(i_row,:) = a_mat(i_row,:) / a_mat(i_row,end);
                end
                if i_metric == 1, scores.CONV_normdbps     = a_mat;   end
                if i_metric == 2, scores.CONV_probperframe = a_mat;   end
            end
            
            glm_conv2{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype} = scores;

        end
    end
end

glm_conv = glm_conv2;
eval(sprintf('save %s/glm_conv.mat glm_conv', glmconvdir))



%}