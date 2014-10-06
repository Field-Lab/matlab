% Dirty Code showing scatter plots of normalized GLM performances
clear; close all; clear all; clc
% Use this for the convergence check!
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS
baseoutput_dir = '/Volumes/Analysis/nora/NSEM/CrossStim_Performance';

% SETUP cells and experiments, the TYPE of GLM (GLMType) 
BD = NSEM_BaseDirectories;
exptests = [1 2 3 4];
cellselectiontype = 'shortlist';%cellselectiontype = 'debug';
GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
%GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
%GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; GLMType.fixedSPlength = 13;  GLMType.fixedSP_nullpoint = 'mean'; 
GLMType.nullpoint = 'mean'; 
%GLMType.fit_type = 'NSEM'; 
GLMType.map_type = 'mapPRJ';
GLMType.debug = false; 
GLMType.specialchange = false;

%GLMType.stimfilter_mode = 'rk1';
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
%GLMType.input_pt_nonlinearity      = true;
%GLMType.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
GLMType.CONVEX = false;
%GLMType.DoubleOpt = true;
%{
GLMType.stimfilter_mode = 'rk1';
GLMType.specialchange = true;
GLMType.specialchange_name = 'ROIlength_9';
GLMType.CONVEX = false;
%}

GLMType.TonicDrive = true;
GLMTYpe.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = true;
GLMType.fixed_spatialfilter = true;
GLMType.func_sname = 'glmwrap24_CP';
GLMType.fullmfilename =mfilename('fullpath'); 
i_exp = 1; i_cell = 1;

GLMType.fitname  = GLM_fitname(GLMType);   
troubleshoot.doit    = true;
troubleshoot.plotdir = '/Volumes/Analysis/nora/NSEM/troubleshootingplot'
%troubleshoot.plotdir = BD.GLM_troubleshootplots  % temporarily change to
troubleshoot.name    = 'singleopt';



agg_perf.byexpnm = cell(4,2)
agg_perf.fitname = GLMType.fitname;
agg_perf.fullGLMType = GLMType;
agg_perf.note1 = 'rows are experiments, (2012-08-09-3,2012-09-27-3,2013-08-19-6,2013-10-10-0)';
agg_perf.note1 = 'columns are fits, (WN, NSEM)';

for i_type = 1:2
    if i_type == 1, GLMType.fit_type = 'WN';  end
    if i_type == 2, GLMType.fit_type = 'NSEM'; end
for i_exp = exptests
    %% 
    expnumber = i_exp;
    [exp_nm,cells,expname]  = cell_list( expnumber, cellselectiontype);
    cells
    agg_perf.byexpnm{i_exp,i_type}.exp_nm = expname;
    agg_perf.byexpnm{i_exp,i_type}.fit_type = GLMType.fit_type; 
    agg_perf.byexpnm{i_exp,i_type}.cells  = cells; 
    agg_perf.byexpnm{i_exp,i_type}.ONP    = zeros(1,length(cells));
    agg_perf.byexpnm{i_exp,i_type}.OFFP   = zeros(1,length(cells));
    agg_perf.byexpnm{i_exp,i_type}.WNperf_vec      = zeros(1,length(cells));  
    agg_perf.byexpnm{i_exp,i_type}.NSEMperf_vec    = zeros(1,length(cells));
    
    agg_perf.byexpnm{i_exp,i_type}.WN_rastsim   = cell(1,length(cells));
    agg_perf.byexpnm{i_exp,i_type}.NSEM_rastsim = cell(1,length(cells));
    
    
    
    [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
    if strcmp(GLMType.fit_type, 'WN')
        SPars_WN = StimulusPars.slv;
        [StimPars ] = Directories_Params_v23(exp_nm, 'NSEM', GLMType.map_type);
        SPars_NSEM  = StimPars.slv;
    elseif strcmp(GLMType.fit_type, 'NSEM')
        SPars_NSEM = StimulusPars.slv;
        [StimPars ] = Directories_Params_v23(exp_nm, 'WN', GLMType.map_type);
        SPars_WN  = StimPars.slv;
    end
        
     %%%%  Shorten Block count if using Debug
    if GLMType.debug
        StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:10);
    end
    clear boolean_debug map_type fit_type shead_cellID expname 
    
    %%%%%% Name and Create a Save Directory %%%%%%%%%%%
        
    inputs.exp_nm    = exp_nm; 
    inputs.map_type  = GLMType.map_type; 
    inputs.stim_type = GLMType.fit_type;
    inputs.fitname   = GLMType.fitname;
  
    d_save = NSEM_secondaryDirectories('savedir_GLMfit', inputs);  clear inputs; 
    display(sprintf('Full Model Fit Parameters are:  %s', GLMType.fitname));  
    display(sprintf('Save Directory :  %s', d_save));
    if ~exist(d_save), mkdir(d_save); end
    GLMType.d_save = d_save; 
    
    %% Load Movie and Concatenate the Fitting Section
    clear Main_SolPars Other_SolParss
    %%% Load Stimulus   -- insert more frame cutting here!    
    %[blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'fitmovie');
    
    [testmovie_WN]   = loadmoviematfile(exp_nm , 'WN', GLMType.cone_model,'testmovie');
    [testmovie_NSEM] = loadmoviematfile(exp_nm , 'NSEM', GLMType.cone_model,'testmovie');
    clear origmatfile
    clear blockedmoviecell blockstartframe fitblocks fitframesperblock framenums
    %}
 
    %% Load Cell Specific Elements   Spikes and STA
    inputs.exp_nm       = exp_nm; 
    inputs.map_type     = GLMType.map_type; 
    DirPars.WN_STAdir   = NSEM_secondaryDirectories('WN_STA', inputs); 

    inputs.stim_type    = 'WN';
    DirPars.organizedspikesdir_WN = NSEM_secondaryDirectories('organizedspikes_dir', inputs);
    
    inputs.stim_type    = 'NSEM';
    DirPars.organizedspikesdir_NSEM = NSEM_secondaryDirectories('organizedspikes_dir', inputs); 
    
    clear inputs
    

    for i_cell = 1:length(cells)
        clear glm_cellstruct
        
        cid = cells{i_cell};
        [celltype , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types); 
	cell_savename
        
        if strcmp(celltype, 'ON-Parasol')
            agg_perf.byexpnm{i_exp,i_type}.ONP(i_cell) = 1;
        elseif strcmp(celltype, 'OFF-Parasol')
            agg_perf.byexpnm{i_exp,i_type}.OFFP(i_cell) = 1;
        else
            error('messed up celltype naming');
        end
        
        outputdir = sprintf('%s/%s/4sec_%s', baseoutput_dir, GLMType.fitname,exp_nm);  
        if ~exist(outputdir, 'dir'), mkdir(outputdir); end
       
        %eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir_WN, cell_savename));
        %organizedspikes_WN = organizedspikes; clear organizedspikes
        eval(sprintf('load %s/%s.mat', d_save, cell_savename));
        
        
        eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir_WN, cell_savename));
        organizedspikes_WN = organizedspikes; clear organizedspikes
        eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir_NSEM, cell_savename));
        organizedspikes_NSEM = organizedspikes; clear organizedspikes;
        % eval(sprintf('load %s/%s.mat', d_save, cell_savename));

	if GLMType.CouplingFilters
            n_couplings=length(fittedGLM.cellinfo.pairs); % number of cells to couple to
            % loading the neighboring spikes to neighborspikes.home
            for j=1:n_couplings
                    % [~ , fittedGLM.cellinfo.pair_savename{j}, ~]  = findcelltype(fittedGLM.cellinfo.pairs(j), datarun_mas.cell_types);

                    	spikesfilename=[DirPars.organizedspikesdir_WN '/organizedspikes_' fittedGLM.cellinfo.pair_savename{j} '.mat'];
			if numel(dir(spikesfilename))
                    		load(spikesfilename,'organizedspikes')
                    		neighbor_organizedspikes_WN{j}=organizedspikes;
			else 
				neighbor_organizedspikes_WN{j}=0;
			end

                    	spikesfilename=[DirPars.organizedspikesdir_NSEM '/organizedspikes_' fittedGLM.cellinfo.pair_savename{j} '.mat'];
			if numel(dir(spikesfilename))
                    		load(spikesfilename,'organizedspikes')
                    		neighbor_organizedspikes_NSEM{j}=organizedspikes;
			else 
				neighbor_organizedspikes_NSEM{j}=0;
			end
            end
        else
                % if there's no coupling, just set this to zero
                neighborspikes=0;neighbor_organizedspikes=0;
        end
        
        fittedGLM0 = fittedGLM;
        fittedGLM = rmfield(fittedGLM,'xvalperformance')

        
        xvalperformance_WN   = eval_xvalperformance_NEW_CP(fittedGLM, SPars_WN, organizedspikes_WN, neighbor_organizedspikes_WN, testmovie_WN)
        xvalperformance_NSEM = eval_xvalperformance_NEW_CP(fittedGLM, SPars_NSEM, organizedspikes_NSEM, neighbor_organizedspikes_NSEM, testmovie_NSEM)
        
        agg_perf.byexpnm{i_exp,i_type}.WNperf_vec(i_cell)   = xvalperformance_WN.glm_normedbits;  
        agg_perf.byexpnm{i_exp,i_type}.NSEMperf_vec(i_cell) = xvalperformance_NSEM.glm_normedbits;
        
        
        agg_perf.byexpnm{i_exp,i_type}.WN_rastsim{i_cell}    = xvalperformance_WN;
        agg_perf.byexpnm{i_exp,i_type}.NSEM_raststim{i_cell} = xvalperformance_NSEM ;
    end
    
end



end
save('agg_perf.mat','agg_perf')



clf;
xlim([0,1]);
ylim([0,1]); hold on;
plot(linspace(0,100,100), linspace(0,100,100),'k' );
set(gca, 'fontsize', 14);
%set(gca,'xtick',[0:20:60]); set(gca,'ytick',[0:20:60]); 
title('GLM Performance out of 100')
MS = 26;
for i_exp = 1:4
    WN_vec   = agg_perf.byexpnm{i_exp,1}.WNperf_vec
    NSEM_vec = agg_perf.byexpnm{i_exp,2}.NSEMperf_vec
    
    a = find(NSEM_vec<0); NSEM_vec(a) = 0; 
    b = find(  WN_vec<0);   WN_vec(b) = 0; 
    
    if i_exp == 1, plot(WN_vec, NSEM_vec,'r.','markersize',MS ); end
    if i_exp == 2, plot(WN_vec, NSEM_vec,'g.','markersize',MS ); end
    if i_exp == 3, plot(WN_vec, NSEM_vec,'b.','markersize',MS ); end
    if i_exp == 4, plot(WN_vec, NSEM_vec,'c.','markersize',MS ); end 
end

xlabel('NSEM-Fit : WN-Test')
ylabel('NSEM-Fit : NSEM-Test')
    
orient landscape
eval(sprintf('print -dpdf NSEMfitsdontcross.pdf'));
%eval(sprintf('print -dpdf WN_outperforms_NSEM.pdf'));










