%%
% AKHEITMAN
% 2015-02-06 Compute Full objective function

% Gives us expTrain_Conv in the convergence analysis data folder

clear; close all; clc
BD = NSEM_BaseDirectories;
basefolder   = sprintf('%s/GLM_Convergence_Analysis/fixedSP_rk1_linear_MU_PS_noCP_timekernelCONEMODEL/Analysis_Plots', BD.NSEM_home);
shortlist = false;
if exist('shortlist','var') && shortlist
    eval(sprintf('load %s/Raw_Conv_shortlist.mat', basefolder));
else
    eval(sprintf('load %s/Raw_Conv.mat', basefolder)); 
end
pct_change = [5:10:95]; pct_change = [pct_change, 100];
celltypes = [1 2]
exps = [1 2 3 4]
i_exp = 1; i_cell = 1; i_celltype = 1; i_pct = 1;  RCD = Raw_Conv_Data;

debug = false;
%%

% CYCLE THROUGH CELLS 
for i_exp = exps
    %%
    exp_nm    = RCD{i_exp}.exp_nm;
    plotdir = sprintf('%s/ConvPlots/%s',basefolder,exp_nm);
    if ~exist(plotdir), mkdir(plotdir); end
    
    % LOAD DIRECTORIES OF STA AND ORGANIZED SPIKES
    
    GLMType   = RCD{i_exp}.GLMType; 
    fitname   = GLM_fitname(GLMType);
    secondDir.exp_nm    = exp_nm; 
	secondDir.map_type  = GLMType.map_type; 
	secondDir.stim_type = 'WN';
	Dirs.WN_organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir); 
	secondDir.stim_type = 'NSEM';
	Dirs.NSEM_organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir); 
    Dirs.WN_STAdir               = NSEM_secondaryDirectories('WN_STA', secondDir);
    clear secondDir
    
    % LOAD STIMULUS PARS AND A DATARUN FOR MASTER
    [StimulusPars,~,~,datarun_mas] = Directories_Params_v23(exp_nm, 'WN', 'mapPRJ');
    SPars.mas = StimulusPars.master;
    SPars.WN  = StimulusPars.slv;
    [StimulusPars] = Directories_Params_v23(exp_nm, 'NSEM', 'mapPRJ');
    SPars.NSEM = StimulusPars.slv;
    clear StimulusPars
    if exist('debug','var') && debug
        SPars.WN.FitBlocks = SPars.WN.FitBlocks(end-1); 
        SPars.NSEM.FitBlocks = SPars.NSEM.FitBlocks(end-1); 
    end
    % HACK FIX
    SPars.NSEM.computedtstim = SPars.WN.computedtstim;
    
    
    % LOAD MOVIES
    [blockedmoviecell, inputstats.WN] = loadmoviematfile(exp_nm , 'WN', GLMType.cone_model,'fitmovie');
    concat_fitmovie.WN      = concat_fitmovie_fromblockedcell(blockedmoviecell, SPars.WN);
    clear blockedmoviecell
    [blockedmoviecell, inputstats.NSEM ] = loadmoviematfile(exp_nm , 'NSEM', GLMType.cone_model,'fitmovie');
    concat_fitmovie.NSEM      = concat_fitmovie_fromblockedcell(blockedmoviecell, SPars.NSEM);
    clear blockedmoviecell
    
    
    %%
    for i_celltype = celltypes
        if i_celltype == 1; cellgroup = RCD{i_exp}.ONP;  celltype = 'ONPar'; end
        if i_celltype == 2; cellgroup = RCD{i_exp}.OFFP; celltype = 'OFFPar'; end
        test_train.WN  = cell(length(cellgroup),1);
        test_train.NSEM  = cell(length(cellgroup),1); 
        
       
        
        %%
        for i_cell = 1:length(cellgroup)
            %%
            
            
            % LOAD SPIKE TIMES
            cid = cellgroup(i_cell); 
            cell_savename = sprintf('%s_%d', celltype,cid);
            display(sprintf('Working on %s: %s', exp_nm, cell_savename));
            
            eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.WN_organizedspikesdir, cell_savename));
            spikesconcat.WN.home = concat_fitspikes_fromorganizedspikes(organizedspikes.block, SPars.WN); 
            eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.NSEM_organizedspikesdir, cell_savename));
            spikesconcat.NSEM.home = concat_fitspikes_fromorganizedspikes(organizedspikes.block, SPars.NSEM);
            
            % LOAD STA AND CENTER COORDINATES 
            eval(sprintf('load %s/STAandROI_%s.mat STAandROI', Dirs.WN_STAdir, cell_savename));
            master_idx         = find(datarun_mas.cell_ids == cid);
            stafit_centercoord = ( datarun_mas.vision.sta_fits{master_idx}.mean );
            stafit_sd          = ( datarun_mas.vision.sta_fits{master_idx}.sd   );
            slvdim.height      = SPars.WN.height; 
            slvdim.width       = SPars.WN.width; 
            [center_coord,sd]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, SPars.mas, slvdim);
            cellinfo.WN.WN_STA = STAandROI.STA;
            cellinfo.WN.slave_centercoord = center_coord;
            slvdim.height      = SPars.NSEM.height; 
            slvdim.width       = SPars.NSEM.width; 
            [center_coord,sd]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, SPars.mas, slvdim);
            cellinfo.NSEM.WN_STA = STAandROI.STA;
            cellinfo.NSEM.slave_centercoord = center_coord;
            cellinfo.NSEM.computedtstim = SPars.NSEM.computedtstim;
            cellinfo.WN.computedtstim = SPars.WN.computedtstim;
            clear center_coord sd stafit_centercoord stafit_sd slvdim
            
            % LOAD FITTED PARAMETERS FROM RAW_CONV_DATA
            if i_celltype == 1
                P_opt.WN   = RCD{i_exp}.CONV_WN_ONP{i_cell}.fit_p;
                P_opt.NSEM = RCD{i_exp}.CONV_NSEM_ONP{i_cell}.fit_p;
                
                XVAL_BPS.WN   = RCD{i_exp}.CONV_WN_ONP{i_cell}.xval_bps;
                XVAL_BPS.NSEM = RCD{i_exp}.CONV_NSEM_ONP{i_cell}.xval_bps;
            elseif i_celltype == 2
                P_opt.WN      = RCD{i_exp}.CONV_WN_OFFP{i_cell}.fit_p;
                P_opt.NSEM    = RCD{i_exp}.CONV_NSEM_OFFP{i_cell}.fit_p;
                
                XVAL_BPS.WN   = RCD{i_exp}.CONV_WN_OFFP{i_cell}.xval_bps;
                XVAL_BPS.NSEM = RCD{i_exp}.CONV_NSEM_OFFP{i_cell}.xval_bps;
            end
            
            % COPMUTE
            [objval_WN, logprob_WN, bps_WN]     = conv_obj_fulltrain_core(GLMType,P_opt.WN, spikesconcat.WN, concat_fitmovie.WN, inputstats.WN, cellinfo.WN);
            [objval_NSEM, logprob_NSEM, bps_NSEM] = conv_obj_fulltrain_core(GLMType,P_opt.NSEM, spikesconcat.NSEM, concat_fitmovie.NSEM, inputstats.NSEM, cellinfo.NSEM);    
            test_train.WN{i_cell}.objval   = objval_WN;
            test_train.NSEM{i_cell}.objval = objval_NSEM;
            test_train.WN{i_cell}.logprob   = logprob_WN;
            test_train.NSEM{i_cell}.logprob = logprob_NSEM;
            test_train.WN{i_cell}.bps   = bps_WN;
            test_train.NSEM{i_cell}.bps = bps_NSEM;

            % PLOT
            clf;
            subplot(3,1,1); 
            set(gca, 'fontsize', 10); axis off
            c = 0;
            c=c+1; text(-.1, 1-0.1*c,sprintf('Exp: %s, Cell %s', exp_nm, cell_savename),'interpreter','none');
            c=c+1; text(-.1, 1-0.1*c,sprintf('Fit Type: %s', fitname),'interpreter','none');
            c=c+1; text(-.1, 1-0.1*c,sprintf('Plot Date %s',datestr(clock)));
            c=c+1; text(-.1, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none' );
            panels = 3;
            for i_fit = 1:2
                if i_fit == 1, minutes = RCD{i_exp}.WN_fit_minutes;   crossval = XVAL_BPS.WN;   traindata =  test_train.WN{i_cell}; fit = 'WN';    string_color = 'b'; end
                if i_fit == 2, minutes = RCD{i_exp}.NSEM_fit_minutes; crossval = XVAL_BPS.NSEM; traindata =  test_train.NSEM{i_cell}; fit = 'NSEM'; string_color = 'r'; end
                for i_panel = 1:panels
                    if i_panel == 1, measure = 'XVAL BPS'; xvals = minutes; yvals = crossval; end
                    if i_panel == 2, measure = 'TRAIN: BPS Trunc0'; xvals = minutes; yvals = traindata.bps; yvals(find(yvals<=0)) =0; end
                    if i_panel == 3, measure = 'TRAIN: -OBJVAL'; xvals = minutes; yvals = -traindata.objval; end
                    subplot(3,panels, (i_panel+(i_fit)*panels) ); hold on;
                    set(gca,'fontsize',8); title(sprintf('%s: %s', fit,measure),'interpreter','none');
                    plot(xvals,yvals,string_color);
                    plot(xvals,yvals,'k.'); 
                    hold off
                end
            end
            orient landscape
            eval(sprintf('print -dpdf %s/CONVPlot_%s.pdf', plotdir, cell_savename))            
        end
        
        if i_celltype == 1, RCD{i_exp}.test_train_ONP  = test_train; end
        if i_celltype == 2, RCD{i_exp}.test_train_OFFP = test_train; end
            
    end
    
    expTrain_Conv = RCD{i_exp};
    if exist('shortlist','var') && shortlist
        eval(sprintf('save %s/expTrain_Conv_%d_shortlist.mat expTrain_Conv GLMType', basefolder,i_exp))
    else
        eval(sprintf('save %s/expTrain_Conv_%d.mat expTrain_Conv GLMType', basefolder,i_exp))
    end
    clear expTrain_Conv
    
end

FullTrain_Conv = RCD;
if exist('shortlist','var') && shortlist
    eval(sprintf('save %s/Train_Conv_shortlist.mat FullTrain_Conv GLMType', basefolder))
else
    eval(sprintf('save %s/Train_Conv.mat FullTrain_Conv GLMType', basefolder))
end
%}