function out = compareFit(dat,fit,sim,plotting)

%-------------------------------------------
%
% out = compareFit(fit,sim)
%
% compare a fitted LN cascade model to
% ground truth (parameters used to simulate responses)
%
% fit must contain the fields:
% fit.S -- number of subunits
% fit.I_sc -- binary assignments from stimulus to subunits, S x C
% fit.A_sc -- weights from stimulus to subunits, S x C
% fit.J_s -- binary assignments from subunits to cell, 1 x S
% fit.B_s -- weights from subunits to cell, 1 x S
%
% sim must contain the fields:
% sim.S -- number of subunits
% sim.I_sc -- binary assignments from stimulus to subunits, S x C
% sim.A_sc -- weights from stimulus to subunits, S x C
% sim.J_s -- binary assignments from subunits to cell, 1 x S
% sim.B_s -- weights from subunits to cell, 1 x S
%
% optional input:
% plotting -- set to 0 to skip plotting (and just
%             return the error) (default = 1)
% 
% this function evaluates the error for all parameters
% and plots the two against one another
%
% freeman, 3-21-2012
%-------------------------------------------

%%
if ~exist('plotting','var') || isempty(plotting)
    plotting = 1;
end

figure; 
set(gcf,'Position',[800 76 600 280]);
set(gcf,'Color','white');

%%
for is=1:fit.S
    inds = find(fit.I_sc(is,:));
    vec1 = sim.A_sc(is,inds);
    vec2 = fit.A_sc(is,inds);
    ang_error_s(is) = subspace(vec1',vec2');
end
nonSingInds = sum(fit.I_sc,2) ~= 1;
out.ang_error_A_sc = mean(ang_error_s(nonSingInds));
out.ang_error_B_s = subspace(fit.B_s',sim.B_s');
%%
%inds = find(fit.I_sc);
%out.mse_A_sc = mean(vector((sim.A_sc(inds)-fit.A_sc(inds)).^2));
%out.mse_B_s = mean((sim.B_s-fit.B_s).^2);

%%
if plotting
%%
subplot(2,4,1);
compVals(sim.B_s(:),fit.B_s(:));
box off
set(gca,'TickDir','out');
xlabel('True'); ylabel('Fitted'); title('Subunit -> RGC weights');


subplot(2,4,2);
nonSingInds = sum(fit.I_sc,2) ~= 1;
tmp1 = sim.A_sc(nonSingInds,:);
tmp2 = fit.A_sc(nonSingInds,:);
compVals(tmp1(find(fit.I_sc(nonSingInds,:))),tmp2(find(fit.I_sc(nonSingInds,:))));
box off
set(gca,'TickDir','out');
xlabel('True'); ylabel('Fitted'); title('Cone -> subunit weights');
%%
subplot(2,4,5);
hold on
plot(fit.B_s,'k.','MarkerSize',20);
plot(sim.B_s,'k-');
box off;
set(gca,'TickDir','out');
axis square
xlim([0 fit.S+1]);
xlabel('Subunit number');
set(gca,'XTick',[1:fit.S]);
label_s{1} = 1;
label_s{fit.S} = fit.S;
set(gca,'XTickLabel',label_s);

nonSingInds = sum(fit.I_sc,2) ~= 1;
nNonSing = sum(nonSingInds);
clrs = zeros(fit.S,3);
clrs(nonSingInds,:) = cmap_hsv_dark(nNonSing);

subplot(2,4,6);
hold on
for is=1:fit.S
    inds = fit.I_sc(is,:)==0;
    tmpfit = fit.A_sc(is,:);
    tmpfit(inds) = NaN;
    tmpsim = sim.A_sc(is,:);
    tmpsim(inds) = NaN;
    if sum(~inds) > 1
    plot(tmpfit,'.','MarkerSize',20,'Color',clrs(is,:));
    plot(sim.A_sc(is,:),'-','Color',clrs(is,:));
    end
end
box off;
set(gca,'TickDir','out');
axis square
xlim([0 fit.C+1]);
ylim([0 max(vector(fit.A_sc(nonSingInds,:)))+0.1]);
xlabel('Cone number');
set(gca,'XTick',[1:fit.C]);
label_c{1} = 1;
label_c{fit.C} = fit.C;
set(gca,'XTickLabel',label_c);
%%
subplot(2,4,3);
imagesc(sim.A_sc);
colormap(gray);
axis off
title('True weights');

subplot(2,4,7);
imagesc(fit.A_sc);
colormap(gray);
axis off
title('Fitted weights');

%%
subplot(2,4,4);
plotNonLin(sim.f);
title('True subunit nonlinearity');

subplot(2,4,8);
plotNonLin(fit.f);
title('Estimated subunit nonlinearity');
end
