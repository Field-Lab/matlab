%{
clear; close all; clc
celltypes    = [1 2]
exps         = [1 2 3 4]
stimtypes    = [1 2]
cell_subset          = 'glmconv_4pct';
changes_cell{1}.type = 'cone_model';
changes_cell{1}.name = 'rieke_linear'
%changes_cell{2}.type = 'input_pt_nonlinearity';
%changes_cell{2}.name = 'piece_linear_aboutmean';
hack_fitintocrossval(exps,stimtypes,celltypes,changes_cell,cell_subset)
%}
function hack_fitintocrossval(exps,stimtypes,celltypes,changes_cell,cell_subset, runoptions)
% DICTATE WHICH EXPERIMENTS AND CELLTYPES TO USE

if exist('runoptions','var')
    if isfield(runoptions,'debug') && runoptions.debug
        debug = true;
    end
end





% SET DIRECTORIES / USING GLM becaause it has the raster.
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection)); 
if exist('changes_cell', 'var')
    GLMType = GLM_settings('default',changes_cell);
else
    GLMType = GLM_settings('default');
end
GLMType.fitname    = GLM_fitname(GLMType)


pars.vkspd_msec         =  50; 
pars.fracvar_smoothbins = 12;


% Victor Bin 
i_exp = 1; i_celltype = 1;  i_cell = 1; i_stimtype = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOP AND COMPUTE
    
%%% REORGANIZE INTO A MORE APPROACHABLE STRUCTURE
for i_exp = exps
    for i_stimtype = stimtypes
        for i_celltype = celltypes
            % Housekeeping
            exp_nm  = allcells{i_exp}.exp_nm;
            expname = allcells{i_exp}.expname;
            if i_stimtype == 1, stimtype = 'WN';   end
            if i_stimtype == 2, stimtype = 'NSEM'; end
            secondDir.exp_nm    = exp_nm;
            secondDir.stim_type = stimtype;
            secondDir.map_type  = 'mapPRJ';
            secondDir.fitname   = GLMType.fitname;
            glmfitdir   = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
            crossdir = sprintf('%s/CrossVal_RawMetrics/', glmfitdir);
            if ~exist(crossdir,'dir'), mkdir(crossdir); end
            map_type= 'mapPRJ';
            if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end

            if strcmp(cell_subset,'all')
                candidate_cells = [allcells{i_exp}.ONP allcells{i_exp}.OFFP]
            elseif strcmp(cell_subset,'shortlist'); 
                [~,candidate_cells,~]  = cell_list(i_exp, 'shortlist'); 
                candidate_cells = cell2mat(candidate_cells) ; 
            elseif strcmp(cell_subset,'glmconv_4pct')
                conv_column = 2; 
                conv_index_ON = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column));
                conv_index_OFF = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column));
                candidate_cells = [allcells{i_exp}.ONP(conv_index_ON) allcells{i_exp}.OFFP(conv_index_OFF)];
            end
            cellgroup = intersect(candidate_cells, cellgroup);
            %cellgroup = fliplr(cellgroup)

            if exist('debug','var') && debug
                cellgroup = cellgroup(1:2);
            end
            %%
            for i_cell = 1:length(cellgroup)
                tStart = tic;
                clear fittedGLM cell_savename
                cid = cellgroup(i_cell); cell_savename = sprintf('%s_%d', celltype,cid);
                
                
                eval(sprintf('load %s/%s.mat', glmfitdir, cell_savename))               
                eval(sprintf('load %s/%s.mat crossval_rawmetrics', crossdir, cell_savename))

                crossval_rawmetrics.fit_objval = fittedGLM.rawfit.objective_val;
                crossval_rawmetrics.note_fitobjval = sprintf('Not CrossValidated.. just a reference fit point, the objective value');
                
                eval(sprintf('save %s/%s.mat crossval_rawmetrics', crossdir, cell_savename))
                duration = toc(tStart);
                display(sprintf('Finished: %s %s %s: Computation Time in Minutes', expname,cell_savename, stimtype,duration/60));
                
            end
           
            
         
        end

    end
end

end


