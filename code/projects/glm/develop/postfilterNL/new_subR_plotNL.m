function new_subR_plotNL(fittedNL_crossprep, plotparams)
% subRoutine form done 2015-06-18 AKHeitman
% PreSet plotting limits, clean colors etc.
% Main Call: only reference to model_comparison structure 
% scores1 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model1;
% scores2 = model_comparison.byexpnm{i_exp}.scores_bycelltype{i_celltype}.finalscores_model2;
% plot_entries

clf
%{
subplot(3,1,1); axis off
delta = .15;
c=0;   text(-.1, 1-delta*c,sprintf('PURPOSE: %s', plotparams.purpose ));
c=c+1; text(-.1, 1-delta*c,sprintf('Metric: %s,  Comparison Title: %s',...
    model_comparison.metric, model_comparison.comparison_name ),'interpreter','none');
c=c+1; text(-.1, 1-delta*c,sprintf('X-axis: %s Fit: %s',...
    model_comparison.models{1}.fit_type,model_comparison.models{1}.fitname),'interpreter','none');
c=c+1; text(-.1, 1-delta*c,sprintf('Y-axis: %s Fit: %s',...
    model_comparison.models{2}.fit_type,model_comparison.models{2}.fitname),'interpreter','none');
c=c+1; text(-.1, 1-delta*c,sprintf('Plot Date %s',datestr(clock)));
c=c+1; text(-.1, 1-delta*c,sprintf('Mfile: %s', mfilename('fullpath')),'interpreter','none' );
%}    


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


% ON only

LW_MU_fat  = 4;
LW_MU_thin = 1;
LW_STD_fat = 1.5;
LW_STD_thin = .5;

clf
% ON Parasols
subplot(3,2,3); clear s_mu s_std SCORES
hold on; set(gca, 'fontsize', 10);
SCORES = ON_scores; title('On Parasols')

axis square
xlabel('Normed Filter Output');
ylabel('Hertz');
s_mu  = mean(SCORES);
s_std = std(SCORES);
xlim([plotparams.minval, plotparams.maxval]);
plot(filter_score,s_mu,'b','linewidth',LW_MU_fat);
plot(filter_score,s_mu,'w','linewidth',LW_MU_thin);
plot(filter_score,s_mu+s_std,'b','linewidth',LW_STD_fat);
plot(filter_score,s_mu+s_std,'w','linewidth',LW_STD_thin);
plot(filter_score,s_mu-s_std,'b','linewidth',LW_STD_fat);
plot(filter_score,s_mu-s_std,'w','linewidth',LW_STD_thin);


% OFF Parasols
subplot(3,2,4); clear s_mu s_std SCORES
hold on; set(gca, 'fontsize', 10);
SCORES = OFF_scores;title('OFF Parasols')

axis square
xlabel('Normed Filter Output');
ylabel('Hertz');
s_mu  = mean(SCORES);
s_std = std(SCORES);
xlim([plotparams.minval, plotparams.maxval]);
plot(filter_score,s_mu,'b','linewidth',LW_MU_fat);
plot(filter_score,s_mu,'k','linewidth',LW_MU_thin);
plot(filter_score,s_mu+s_std,'b','linewidth',LW_STD_fat);
plot(filter_score,s_mu+s_std,'k','linewidth',LW_STD_thin);
plot(filter_score,s_mu-s_std,'b','linewidth',LW_STD_fat);
plot(filter_score,s_mu-s_std,'k','linewidth',LW_STD_thin);


% Combine On and Off
subplot(3,2,5); clear s_mu s_std SCORES
hold on; set(gca, 'fontsize', 10);
SCORES = [ON_scores; OFF_scores]; title('ON OFF combined')

axis square
xlabel('Normed Filter Output');
ylabel('Hertz');
s_mu  = mean(SCORES);
s_std = std(SCORES);
xlim([plotparams.minval, plotparams.maxval]);
plot(filter_score,s_mu,'b','linewidth',LW_MU_fat);
plot(filter_score,s_mu+s_std,'b','linewidth',LW_STD_fat);
plot(filter_score,s_mu-s_std,'b','linewidth',LW_STD_fat);


subplot(3,2,6); clear s_mu s_std SCORES
hold on; set(gca, 'fontsize', 10);
title('ON OFF Seperate')

axis square
xlabel('Normed Filter Output');
ylabel('Hertz');
SCORES = ON_scores;
s_mu  = mean(SCORES);
s_std = std(SCORES);
xlim([plotparams.minval, plotparams.maxval]);
plot(filter_score,s_mu,'b','linewidth',LW_MU_fat);
plot(filter_score,s_mu,'w','linewidth',LW_MU_thin);
%{
plot(filter_score,s_mu+s_std,'b','linewidth',LW_STD_fat);
plot(filter_score,s_mu+s_std,'w','linewidth',LW_STD_thin);
plot(filter_score,s_mu-s_std,'b','linewidth',LW_STD_fat);
plot(filter_score,s_mu-s_std,'w','linewidth',LW_STD_thin);
%}

SCORES = OFF_scores;
s_mu  = mean(SCORES);
s_std = std(SCORES);
xlim([plotparams.minval, plotparams.maxval]);
plot(filter_score,s_mu,'b','linewidth',LW_MU_fat);
plot(filter_score,s_mu,'k','linewidth',LW_MU_thin);

%{
plot(filter_score,s_mu+s_std,'b','linewidth',LW_STD_fat);
plot(filter_score,s_mu+s_std,'k','linewidth',LW_STD_thin);
plot(filter_score,s_mu-s_std,'b','linewidth',LW_STD_fat);
plot(filter_score,s_mu-s_std,'k','linewidth',LW_STD_thin);
%}


orient tall
eval(sprintf('print -dpdf plot_%s.pdf',printname))
end