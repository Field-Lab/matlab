%{

clear; clc; homedir = pwd
exps = [1 2 3 4]; 
cell_subset = 'glmconv_4pct'; 
%baseGLM.settings{1}.type = 'PostSpikeFilter';
%baseGLM.settings{1}.name =  'OFF';
baseGLM.settings  = {};
%baseGLM.settings{1}.type = 'cone_model';
%baseGLM.settings{1}.name = 'rieke_linear'
%baseGLM.settings{2}.type= 'input_pt_nonlinearity';
%baseGLM.settings{2}.name= 'piecelinear_fourpiece_eightlevels';
postfilterNL.type = 'Logistic_fixMU_noPS';

for i_stimtypes = [2 1]
    if i_stimtypes == 1, stimtype = 'WN'; end
    if i_stimtypes == 2, stimtype = 'NSEM'; end
    plot_fittedNL(exps,stimtype,cell_subset,baseGLM.settings,postfilterNL)%runoptions)
end

%}
function plot_fittedNL(exps,stimtype,cell_subset,baseGLM_settings,postfilterNL)
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
Dirs.save    = sprintf('%s/Plots/Fitted_NL/%s/%s', BD.GLM_output_analysis, postfilterNL.type, baseGLM.Type.fitname);
if ~exist(Dirs.save, 'dir'), mkdir(Dirs.save); end


% Compute and Save Means and STDS  (OR LOAD IF ALREADY DONE)
fittedNL_crossprep = cell(4,1);
for i_exp = exps    
    % Load master datarun, bookkeep
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    
    savename    = sprintf('%s_%s_%s',stimtype,cell_subset,exp_nm);
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




% Plot HERE\
plotparams.exps = exps;
plotparams.name_NL = postfilterNL.type;
plotparams.underlying_model = baseGLM.Type.fitname;

expstring = 'exps';
if find(plotparams.exps ==1), expstring = sprintf('%sA',expstring); end
if find(plotparams.exps ==2), expstring = sprintf('%sB',expstring); end
if find(plotparams.exps ==3), expstring = sprintf('%sC',expstring); end
if find(plotparams.exps ==4), expstring = sprintf('%sD',expstring); end
printname.base = sprintf('%s_%s_%s',stimtype,cell_subset,expstring);
printname.separate = sprintf('%s_separate',printname.base);
printname.combined = sprintf('%s_combined',printname.base);


Dirs.print = sprintf('%s/Plots/', Dirs.save);
if ~exist(Dirs.print, 'dir'), mkdir(Dirs.print); end
cd(Dirs.print)
for domainval = [3 5 7 9]
    plotparams.minval = -domainval;
    plotparams.maxval =  domainval;
    
    plotparams.type   = 'combined';
    printname.final = sprintf('%s_Dom%d', printname.combined,domainval);
    subR_plotNL(fittedNL_crossprep, plotparams, printname.final)

    plotparams.type   = 'separate';
    printname.final   = sprintf('%s_Dom%d',printname.separate,domainval);
    subR_plotNL(fittedNL_crossprep, plotparams, printname.final)
end
%

% hack to get plot without the 2013-10-10-0
if strcmp(expstring, 'expsABCD')
    no4.expstring  = 'expsABC';
    no4.printname.base = sprintf('%s_%s_%s',stimtype,cell_subset,no4.expstring);
    no4.printname.separate = sprintf('%s_separate',no4.printname.base);
    no4.printname.combined = sprintf('%s_combined',no4.printname.base);

    for domainval = [3 5 7 9]
        plotparams.minval = -domainval;
        plotparams.maxval =  domainval;

        plotparams.type      = 'combined';
        no4.printname.final  = sprintf('%s_Dom%d', no4.printname.combined,domainval);
        subR_plotNL(fittedNL_crossprep, plotparams, no4.printname.final)

        plotparams.type       = 'separate';
        no4.printname.final   = sprintf('%s_Dom%d',no4.printname.separate,domainval);
        subR_plotNL(fittedNL_crossprep, plotparams, no4.printname.final)
    end
    clear no4
end


plotparams0 = plotparams;
clear plotparams;
for i_exp = 1:length(plotparams0.exps)
    plotparams = plotparams0;
    plotparams.exps = plotparams.exps(i_exp);
    expstring = 'exps';
    if find(plotparams.exps ==1), expstring = sprintf('%sA',expstring); end
    if find(plotparams.exps ==2), expstring = sprintf('%sB',expstring); end
    if find(plotparams.exps ==3), expstring = sprintf('%sC',expstring); end
    if find(plotparams.exps ==4), expstring = sprintf('%sD',expstring); end
    printname.base = sprintf('%s_%s_%s',stimtype,cell_subset,expstring);
    printname.separate = sprintf('%s_separate',printname.base);
    printname.combined = sprintf('%s_combined',printname.base);
    for domainval = [3 5 7 9]
        plotparams.minval = -domainval;
        plotparams.maxval =  domainval;

        plotparams.type   = 'combined';
        printname.final = sprintf('%s_Dom%d', printname.combined,domainval);
        subR_plotNL(fittedNL_crossprep, plotparams, printname.final)

        plotparams.type   = 'separate';
        printname.final   = sprintf('%s_Dom%d',printname.separate,domainval);
        subR_plotNL(fittedNL_crossprep, plotparams, printname.final)
    end
end
cd(Dirs.home)



end

%%
function [output_LOGI] = subR_compute_LOGISTIC(MAX,RATE, Y_INT, input_LOGI)
    OFFSET   = log( (MAX/Y_INT) - 1  ) / RATE;     
    output_LOGI = (MAX ./ (1 + exp(-RATE * (input_LOGI- OFFSET) )));
end
function subR_plotNL(fittedNL_crossprep, plotparams, printname)
% subRoutine form done 2015-06-18 AKHeitman
% PreSet plotting limits, clean colors etc.
% Main Call: only reference to model_comparison structure 
% scores1 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model1;
% scores2 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model2;
% plot_entries



%{
printname = 'yeah';
plotparams.minval = -9;
plotparams.maxval = 9;
plotparams.exps = [1 2 3];
plotparams.type = 'combined';
subR_plotNL(fittedNL_crossprep, plotparams, printname)

%}
clf
%
subplot(3,1,1); axis off
delta = .15;
c=0;   text(-.1, 1-delta*c,sprintf('NONLINEARITY: %s', plotparams.name_NL ),'interpreter','none');
c=c+1; text(-.1, 1-delta*c,sprintf('PreNL GLM fit: %s',plotparams.underlying_model),'interpreter','none');
c=c+1; text(-.1, 1-delta*c,sprintf('X-axis: Output of Linear Stimulus Filter Z-scored individually for each cell'));
c=c+1; text(-.1, 1-delta*c,sprintf('White Center Stipe: On Parasol, Black Center Stripe: Off Parasol'));
c=c+1; text(-.1, 1-delta*c,sprintf('Plot Date %s',datestr(clock)));
c=c+1; text(-.1, 1-delta*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none' );
%}    


main_color = 'b';
highlight_ON = 'w';
highlight_OFF = 'k';

% Collect the Data
ON_scores  = [];
OFF_scores = [];
for i_exp = plotparams.exps
    NL = fittedNL_crossprep{i_exp};
    filter_score = NL.ONP_fits.input;
    new_ON_scores = cell2mat(NL.ONP_fits.output_opt_bycell);
    new_OFF_scores = cell2mat(NL.OFFP_fits.output_opt_bycell);
    
    ON_scores = [ON_scores ; new_ON_scores];
    OFF_scores = [OFF_scores ; new_OFF_scores];
end

domain_ind = intersect( find(filter_score<=plotparams.maxval) , find(filter_score>=plotparams.minval) );
filter_score = filter_score(domain_ind);
ON_scores = ON_scores(:,domain_ind);
OFF_scores = OFF_scores(:,domain_ind);
%filter_score = filter_score(find(filter_score<=plotparams.maxval));
%filter_score = filter_score(find(filter_score>=plotparams.minval));

% ON only
LW_MU_fat  = 8;
LW_MU_thin = 3;
LW_STD_fat = 1;
LW_STD_thin = .5;

y_min = 0;
y_max = max( [max(mean(ON_scores)+std(ON_scores)),  max(mean(OFF_scores)+std(OFF_scores))] );

if strcmp(plotparams.type, 'separate')

    % ON Parasols
    subplot(3,2,[3 5]);
    hold on; set(gca, 'fontsize', 10);
    SCORES = ON_scores; title('On Parasols')
    xlim([plotparams.minval, plotparams.maxval]);
    ylim([y_min,y_max]);
    
    
    axis square
    xlabel('Stim Filter Output (z score)');
    ylabel('Hertz');
    s_mu  = mean(SCORES);
    s_std = std(SCORES);
    xlim([plotparams.minval, plotparams.maxval]);
    plot(filter_score,s_mu,main_color,'linewidth',LW_MU_fat);
    plot(filter_score,s_mu,highlight_ON,'linewidth',LW_MU_thin);
    plot(filter_score,s_mu+s_std,main_color,'linewidth',LW_STD_fat);
    plot(filter_score,s_mu+s_std,highlight_ON,'linewidth',LW_STD_thin);
    plot(filter_score,s_mu-s_std,main_color,'linewidth',LW_STD_fat);
    plot(filter_score,s_mu-s_std,highlight_ON,'linewidth',LW_STD_thin);


    % OFF Parasols
    subplot(3,2,[4 6]); clear s_mu s_std SCORES
    hold on; set(gca, 'fontsize', 10);
    SCORES = OFF_scores;title('OFF Parasols')
    xlim([plotparams.minval, plotparams.maxval]);
    ylim([y_min,y_max]);
    axis square
    xlabel('Stim Filter Output (z score)');
    ylabel('Hertz');
    s_mu  = mean(SCORES);
    s_std = std(SCORES);
    xlim([plotparams.minval, plotparams.maxval]);
    plot(filter_score,s_mu,main_color,'linewidth',LW_MU_fat);
    plot(filter_score,s_mu,highlight_OFF,'linewidth',LW_MU_thin);
    plot(filter_score,s_mu+s_std,main_color,'linewidth',LW_STD_fat);
    plot(filter_score,s_mu+s_std,highlight_OFF,'linewidth',LW_STD_thin);
    plot(filter_score,s_mu-s_std,main_color,'linewidth',LW_STD_fat);
    plot(filter_score,s_mu-s_std,highlight_OFF,'linewidth',LW_STD_thin);
    
    orient landscape
    eval(sprintf('print -dpdf %s.pdf',printname)) 
end

if strcmp(plotparams.type, 'combined')
    % Combine On and Off
    subplot(3,2,[3 5]); clear s_mu s_std SCORES
    hold on; set(gca, 'fontsize', 10);
    SCORES = [ON_scores; OFF_scores]; title('ON OFF combined')
    xlim([plotparams.minval, plotparams.maxval]);
    ylim([y_min,y_max]);

    axis square
    xlabel('Stim Filter Output (z score)');
    ylabel('Hertz');
    s_mu  = mean(SCORES);
    s_std = std(SCORES);
    xlim([plotparams.minval, plotparams.maxval]);
    plot(filter_score,s_mu,main_color,'linewidth',LW_MU_fat);
    plot(filter_score,s_mu+s_std,main_color,'linewidth',LW_STD_fat);
    plot(filter_score,s_mu-s_std,main_color,'linewidth',LW_STD_fat);


    subplot(3,2,[4 6]); clear s_mu s_std SCORES
    hold on; set(gca, 'fontsize', 10);
    title('ON OFF separate')
    xlim([plotparams.minval, plotparams.maxval]);
    ylim([y_min,y_max]);

    axis square
    xlabel('Stim Filter Output (z score)');
    ylabel('Hertz');
    SCORES = ON_scores;
    s_mu  = mean(SCORES);
    s_std = std(SCORES);
    xlim([plotparams.minval, plotparams.maxval]);
    plot(filter_score,s_mu,main_color,'linewidth',LW_MU_fat);
    plot(filter_score,s_mu,highlight_ON,'linewidth',LW_MU_thin);
    
    plot(filter_score,s_mu+s_std,main_color,'linewidth',LW_STD_fat);
    plot(filter_score,s_mu+s_std,highlight_ON,'linewidth',LW_STD_thin);
    plot(filter_score,s_mu-s_std,main_color,'linewidth',LW_STD_fat);
    plot(filter_score,s_mu-s_std,highlight_ON,'linewidth',LW_STD_thin);
    

    SCORES = OFF_scores;
    s_mu  = mean(SCORES);
    s_std = std(SCORES);
    xlim([plotparams.minval, plotparams.maxval]);
    plot(filter_score,s_mu,main_color,'linewidth',LW_MU_fat);
    plot(filter_score,s_mu,highlight_OFF,'linewidth',LW_MU_thin);

    
    plot(filter_score,s_mu+s_std,main_color,'linewidth',LW_STD_fat);
    plot(filter_score,s_mu+s_std,highlight_OFF,'linewidth',LW_STD_thin);
    plot(filter_score,s_mu-s_std,main_color,'linewidth',LW_STD_fat);
    plot(filter_score,s_mu-s_std,highlight_OFF,'linewidth',LW_STD_thin);
    


    orient landscape
    eval(sprintf('print -dpdf %s.pdf',printname))
end

end