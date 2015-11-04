function new_subR_plotNL(fittedNL_crossprep, plotparams, printname)
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
new_subR_plotNL(fittedNL_crossprep, plotparams, printname)

%}
clf
%{
subplot(3,1,1); axis off
delta = .15;
c=0;   text(-.1, 1-delta*c,sprintf('PURPOSE: %s', plotparams.purpose ));
c=c+1; text(-.1, 1-delta*c,sprintf('X-axis: Output of Linear Stimulus Filter
c=c+1; text(-.1, 1-delta*c,sprintf('Y-axis: %s Fit: %s',...
    model_comparison.models{2}.fit_type,model_comparison.models{2}.fitname),'interpreter','none');
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
    eval(sprintf('print -dpdf plot_%s.pdf',printname)) 
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