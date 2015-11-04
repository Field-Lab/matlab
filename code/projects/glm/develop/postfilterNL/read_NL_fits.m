%{

clear; clc
exps = [1 2 3 4]; stimtypes = 'NSEM'; 
cell_subset = 'glmconv_4pct'; 
%baseGLM.settings{1}.type = 'PostSpikeFilter';
%baseGLM.settings{1}.name =  'OFF';
baseGLM.settings = {};
postfilterNL.type        = 'Logistic_fixMU_noPS';
read_NL_fits(exps,stimtypes,cell_subset,baseGLM.settings,postfilterNL)%runoptions)

%}
function read_NL_fits(exps,stimtype,cell_subset,baseGLM_settings,postfilterNL,runoptions)
% started: AKHeitman 2015-06-25

% Load core directories and all eligible cells
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));
% Define structure which uniquely defines GLM to be used 
if exist('baseGLM_settings', 'var')
    baseGLM.settings = baseGLM_settings; clear baseGLM_settings;
    baseGLM.Type = GLM_settings('default',baseGLM.settings);
else
    baseGLM.Type = GLM_settings('default');
end
baseGLM.Type.fitname    = GLM_fitname(baseGLM.Type); 

Dirs.basedir = sprintf('%s/%s/%s', BD.GLM_output_analysis, ...
    baseGLM.Type.fitname, postfilterNL.type);

Dirs.home    = pwd;
Dirs.save    = sprintf('%s/Plots/Fitted_NL/%s', BD.GLM_output_analysis, postfilterNL.type);
if ~exist(Dirs.save, 'dir'), mkdir(Dirs.save); end


% Compute and Save Means and STDS  (OR LOAD IF ALREADY DONE)
fittedNL_crossprep = cell(4,1);
for i_exp = exps    
    % Load master datarun, bookkeep
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    
    savename    = sprintf('%s_%s',cell_subset,exp_nm);
    matfilename = sprintf('%s/%s.mat',Dirs.save,savename);
    
    if exist(matfilename)
        eval(sprintf('load %s/%s.mat fittedNL', Dirs.save, savename));
        fittedNL_crossprep{i_exp} = fittedNL;
    else
        Dirs.loadfits = sprintf('%s/%s_mapPRJ/%s', Dirs.basedir,stimtype,exp_nm);
        fittedNL = allcells{i_exp};
        % Loop through cells 
        for i_celltype = [1 2];            
            if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
            if strcmp(cell_subset,'all')
                candidate_cells = [allcells{i_exp}.ONP allcells{i_exp}.OFFP]
            elseif strcmp(cell_subset,'shortlist') || strcmp(cell_subset, 'debug') 
                [~,candidate_cells,~]  = cell_list(i_exp, cell_subset); 
                candidate_cells = cell2mat(candidate_cells) ; 
            elseif strcmp(cell_subset,'glmconv_4pct')
                eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection));              
                conv_column = 2; 
                conv_index_ON = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column));
                conv_index_OFF = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column));
                candidate_cells = [allcells{i_exp}.ONP(conv_index_ON) allcells{i_exp}.OFFP(conv_index_OFF)];
            end
            cellgroup = intersect(candidate_cells, cellgroup);

            NL.input  = linspace(-10,10,1000);
            NL.output_opt_bycell = cell(length(cellgroup),1); 
            NL.output_standardGLM_bycell = cell(length(cellgroup),1);
            NL.cellgroup = cellgroup;

            for i_cell = 1:length(cellgroup)
                %% Actual Computation
                cid = cellgroup(i_cell); 
                cell_savename = sprintf('%s_%d', celltype,cid);
                display(sprintf('working on %s: %s', expname, cell_savename))            

                if strcmp(postfilterNL.type,'Logistic_fixMU_noPS')
                    eval(sprintf('load %s/%s.mat fittedGLM', Dirs.loadfits, cell_savename))
                    MAX           = fittedGLM.NL_Output.maxrate;
                    RATE          = fittedGLM.NL_Output.slope;
                    Y_INT         = fittedGLM.stimtransform.inputNL.y_int;
                    [output_NL]   = subR_compute_LOGISTIC(MAX,RATE,Y_INT,NL.input);
                    scale_forexp        = fittedGLM.linearfilters.Stimulus_rescale;
                    output_standardGLM = exp(scale_forexp*NL.input)*Y_INT;
                end

                NL.output_opt_bycell{i_cell} = output_NL;
                NL.output_standardGLM_bycell{i_cell} = output_standardGLM;
            end
            NL.opt_mean = mean(cell2mat(NL.output_opt_bycell));
            NL.opt_std  = std(cell2mat(NL.output_opt_bycell));
            NL.standardGLM_mean = mean(cell2mat(NL.output_standardGLM_bycell));
            NL.standardGLM_std  = std(cell2mat(NL.output_standardGLM_bycell));
            
            if i_celltype == 1, fittedNL.ONP_fits = NL; end
            if i_celltype == 2, fittedNL.OFFP_fits = NL; end
        end
        eval(sprintf('save %s/%s.mat fittedNL', Dirs.save, savename));
        fittedNL_crossprep{i_exp} = fittedNL;
    end
end


% Plot HERE
%{
plotparams.exps = exps;
expstring = 'exps';
plotname = sprintf('plots_%s', cell_subset)
if find(plotparams.exps ==1), expstring = sprintf('%sA',expstring); end
if find(plotparams.exps ==2), expstring = sprintf('%sB',expstring); end
if find(plotparams.exps ==3), expstring = sprintf('%sC',expstring); end
if find(plotparams.exps ==4), expstring = sprintf('%sD',expstring); end
printname_base = sprintf('%s_%s',plotname,expstring);
printname_notext = sprintf('%s_NOTEXT',printname_base);
printname_fullplot = sprintf('%s_FULLPLOTS',printname_base);

cd(savedir)
subR_plotcomparison_nolabels(model_comparison,plotparams,printname_notext);
subR_plotcomparison_fullplots(model_comparison,plotparams,printname_fullplot);
cd(homedir)

% hack to get plot without the 2013-10-10-0
if strcmp(expstring, 'expsABCD')
    expstring_no4  = 'expsABC';
    plotparams_no4 = plotparams;
    plotparams_no4.exps = [2 3 1];
    printname_notext   = sprintf('%s_%s_NOTEXT',plotname,expstring_no4);
    printname_fullplot = sprintf('%s_%s_FULLPLOTS',plotname,expstring_no4);
    
    cd(savedir)
    subR_plotcomparison_nolabels(model_comparison,plotparams_no4,printname_notext);
    subR_plotcomparison_fullplots(model_comparison,plotparams_no4,printname_fullplot);
    cd(homedir)
    clear expstring_no4 plotparams_no4
end

subset_params = plotparams;
for i_exp = 1:length(plotparams.exps)
    subset_params.exps = plotparams.exps(i_exp);
    expstring = 'exps';
    if find(subset_params.exps ==1), expstring = sprintf('%sA',expstring); end
    if find(subset_params.exps ==2), expstring = sprintf('%sB',expstring); end
    if find(subset_params.exps ==3), expstring = sprintf('%sC',expstring); end
    if find(subset_params.exps ==4), expstring = sprintf('%sD',expstring); end

    for i_celltype = 1:2
        subset_params.celltypes = i_celltype;
        if i_celltype == 1, celltypestring = 'ONParasols'; end
        if i_celltype == 2, celltypestring = 'OFFParasols'; end
        printname_notext = sprintf('%s_%s_%s_NOTEXT',plotname,expstring,celltypestring);
        printname_fullplot = sprintf('%s_%s_%s_FULLPLOTS',plotname,expstring,celltypestring);
        cd(savedir)
        subR_plotcomparison_nolabels(model_comparison,subset_params,printname_notext);
        subR_plotcomparison_fullplots(model_comparison,subset_params,printname_fullplot);
        cd(homedir)
        
    end
end
%}




end
function [output_LOGI] = subR_compute_LOGISTIC(MAX,RATE, Y_INT, input_LOGI)
    OFFSET   = log( (MAX/Y_INT) - 1  ) / RATE;     
    output_LOGI = (MAX ./ (1 + exp(-RATE * (input_LOGI- OFFSET) )));
end

function subR_plotfittedNL_pop(fittedGLM, fittedGLM_preNL, savedir)

% Cleaned up AKHeitman 2015-06-24
homedir = pwd;
clf;  
printname = sprintf('DiagNLPlot_%s',fittedGLM.cell_savename);
info    = fittedGLM.cellinfo;
GLMType = fittedGLM_preNL.GLMType;

% text box
subplot(3,1,1)
axis off
set(gca, 'fontsize', 12)
obj_NEW = fittedGLM.rawfit.objective_val;
obj_OLD = fittedGLM_preNL.rawfit.objective_val;
optNL_describe  = fittedGLM.NL_Output.note_metric;
optNL_string    = fittedGLM.NL_Output.param_string;

c = 0; offset = 0; delta_c = 1.1;
text(-offset, 1-0.1*c,sprintf('%s: %s %d: %s-Fit (red): POSTNL refit with %s',...
    info.exp_nm, info.celltype,info.cid, GLMType.fit_type, fittedGLM.nonlinearity), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('Red is original GLM, Blue includes Postfilter Nonlinearity'))
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('Objval PctChange %d: from %1.2e to %1.2e',...
   round(100*(obj_NEW-obj_OLD)/obj_OLD),obj_OLD,obj_NEW), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('%s', optNL_describe), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('%s', optNL_string), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('Original Fit: %s',GLMType.fitname), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('NL fit Computated at %s',datestr(clock)), 'interpreter','none')
c = c + delta_c; 
text(-offset, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')), 'interpreter','none' );



% plotting non-linearities
x1 = sort(fittedGLM.stimtransform.normalized_filteroutput_fit); 
y1 = sort(fittedGLM.stimtransform.cif_rawGLM_fit);
y2 = sort(fittedGLM.stimtransform.cif_withNL_fit);

LW = 2;
subplot(3,3,4); set(gca, 'fontsize', 10); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([-4,4]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity - Central Portion')

subplot(3,3,5); set(gca, 'fontsize', 10); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([ max(min(x1),-10),0]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity: Inhibitory Portion')

subplot(3,3,6); set(gca, 'fontsize', 10); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([0,min(10,max(x1))]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity: Excitatory Portion')


% plot rasters
subplot(3,1,3); set (gca, 'fontsize',10)
secs     = 6;
dt = fittedGLM.t_bin;
bins     = 120 * 8 * fittedGLM.bins_per_frame;
rec_rast = fittedGLM.xvalperformance.rasters.recorded(:,1:bins);
glm_rast = fittedGLM_preNL.xvalperformance.rasters.glm_sim(:,1:bins);  
NL_rast  = fittedGLM.xvalperformance.rasters.glm_sim(:,1:bins); 
trials   = size(rec_rast,1);
time     = dt*[1:bins];
xlim([0, ceil(time(end))]);
ylim([0 , 3*trials]); hold on
for i_trial = 1:trials
    rec1 = time(find(rec_rast(i_trial,:)));
    glm1 = time(find(glm_rast(i_trial,:)));
    NL1  = time(find(NL_rast(i_trial,:)));
    % Plot the raster
    plot(rec1, i_trial, 'k.')


    yshift = i_trial;
    if length(glm1) < 4*length(rec1) 
        if length(glm1) > 0
            plot(glm1, yshift + trials, 'r.')
        end
    end
    if length(NL1) < 4*length(rec1) 
        if length(NL1) > 0
            plot(NL1, yshift + 2*trials, 'b.')
        end
    end
end
xlabel('seconds'); ylabel('trials')

cd(savedir)
orient landscape
eval(sprintf('print -dpdf %s.pdf',printname))
cd(homedir)
end       


