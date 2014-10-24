%% set up the simulation
clear all
%%
nSim = 5;
durs = [1 2 4 8 16 32]*60*12;
%durs = 16*60*12;
%err_A_sc = zeros(nSim,length(durs));
%err_B_s = zeros(nSim,length(durs));
for isim=11:25
for idur = 1:length(durs)
    % stimulus parameters
    stim.type = 'gaussian';
    stim.len = [durs(idur) 10000];
    stim.p = 1;
    
    % weight and assignment parameters
    sim.type = 'parasol-overlap';
    sim.C = 10;
    [sim.I_sc sim.locs_c] = mkSubAssgn(sim.C,sim.type); % cone to subunit assignments
    sim.S = size(sim.I_sc,1);
    sim.C = size(sim.I_sc,2);
    sim.J_s = ones(1,sim.S); % subunit to RGC assignments
    
    load parasol-overlap-sim.mat
    sim.B_s = fit.B_s;%ones(size(fit.B_s));
    sim.A_sc = fit.A_sc;
    sim.f.type = 'rectLinear';
    sim.f.p = [0 1];
    sim.g = fit.g;
    clear fit;
    
    % generate stimuli
    stim.C = sim.C;
    [train test] = getStim(stim); % stimulus matrix
    
    % get the training responses
    Y_st = evalNonLin((sim.I_sc.*sim.A_sc)*train.X_ct,sim.f);
    Z_t = evalNonLin((sim.J_s.*sim.B_s)*Y_st,sim.g);
    train.R_t = poissrnd(Z_t);
    
    % get the test responses
    Y_st = evalNonLin((sim.I_sc.*sim.A_sc)*test.X_ct,sim.f);
    Z_t = evalNonLin((sim.J_s.*sim.B_s)*Y_st,sim.g);
    test.R_t = poissrnd(Z_t);
    
    % fitting options
    fit.error = 'loglik';
    fit.constraints.B_s = 'none';
    fit.constraints.A_sc = 'L1';
    fit.S = size(sim.I_sc,1);
    fit.C = size(sim.I_sc,2);
    fit.locs_c = sim.locs_c;
    
    %% try to infer both
    % stuff we'll assume
    fit.I_sc = sim.I_sc;
    fit.J_s = ones(1,size(fit.I_sc,1));
    fit.S = size(fit.I_sc,1);
    fit.C = size(fit.I_sc,2);
    
    % stuff we are initializing
    fit.g.type = 'logexp1';
    fit.g.p = [4 1];
    %fit.f.type = 'spline';
    %fit.f.init = 'rectLinear';
    %fit.f.nknots = 8;
    %fit.f.smoothness = 3;
    %fit.f.extrap = 1;
    %fit.f.w = [];
    fit.f.type = 'rectLinear';
    fit.f.p = [0 1];
    fit.B_s = normMat(ones(1,fit.S),'L1',2);
    fit.A_sc = normMat(rand(fit.S,fit.C).*fit.I_sc,'L1',2);
    %fit = fitF(train,fit,1);
    %%
    %plotModel(fit); drawnow;
    nIter = 5;
    for iIter = 1:nIter
        fit = fitB(train,fit,0);
        fit = fitA(train,fit,0);
        %plotModel(fit); drawnow;
    end
    tmp = compareFit(fit,sim,0);
    err_A_sc(isim,idur) = tmp.ang_error_A_sc;
    err_B_s(isim,idur) = tmp.ang_error_B_s;
    err_A_sc
    err_B_s
end
end
%%
err = compareFit(fit,sim);

%%
clf
set(gcf,'Position',[500   400   730   330]);
clrs = cmap_hsv_dark(3);
clrs = [0.2 0.2 0.2; 0.5 0.5 0.5; 0.8 0.8 0.8];
grays = bsxfun(@times,ones(6,3),linspace(0.9,0.1,6)');
subplot(1,2,1);
hold on
for i=1:3
h = plot(durs/(60*12),mean(results(i).err_A_sc),'-');
set(h,'Color',clrs(i,:));
set(h,'MarkerSize',40);
set(h,'MarkerFaceColor','none');
set(h,'MarkerEdgeColor','k');
set(h,'LineWidth',3);
tmp = mean(results(i).err_A_sc);
for j=1:length(tmp)
    h2 = plot(durs(j)/(60*12),tmp(j),'.');
    set(h2,'MarkerSize',40);
    set(h2,'MarkerFaceColor','white');
    set(h2,'MarkerEdgeColor',grays(j,:));
end
end
set(gca,'XScale','log');
xlim([0.8 38]);
ylim([0 1]);
box off;
axis square
set(gca,'TickDir','out');
ylabel('Angular error','FontSize',16,'FontAngle','italic');
xlabel('Experiment duration (min)','FontSize',16,'FontAngle','italic');
set(gca,'YTick',[0:0.2:1]);
set(gca,'YTickLabel',[0:0.2:1]);
set(gca,'FontSize',12);
title('Cone to subunit weights','FontSize',20,'FontAngle','italic');

subplot(1,2,2);
hold on
for i=1:3
h = plot(durs/(60*12),mean(results(i).err_B_s),'-');
set(h,'MarkerSize',40);
set(h,'MarkerFaceColor','none');
set(h,'MarkerEdgeColor','k');
set(h,'LineWidth',3);
set(h,'Color',clrs(i,:));
tmp = mean(results(i).err_B_s);
for j=1:length(tmp)
    h2 = plot(durs(j)/(60*12),tmp(j),'.');
    set(h2,'MarkerSize',40);
    set(h2,'MarkerFaceColor','white');
    set(h2,'MarkerEdgeColor',grays(j,:));
end
end
set(gca,'XScale','log');
xlim([0.8 38]);
ylim([0 1]);
box off;
axis square
set(gca,'TickDir','out');
ylabel('Angular error','FontSize',16,'FontAngle','italic');
xlabel('Experiment duration (min)','FontSize',16,'FontAngle','italic');
set(gca,'YTick',[0:0.2:1]);
set(gca,'YTickLabel',[0:0.2:1]);
set(gca,'FontSize',14);
title('Subunit to RGC weights','FontSize',20,'FontAngle','italic');

