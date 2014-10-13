%{

for pair=1:6

        set(gca, 'fontsize', 10);
        axis off;
plot_rf_fit(datarun_mas, glm_cellinfo.pairs(1:12),'edge',true)
plot_rf_fit(datarun_mas, glm_cellinfo.pairs(pair), 'fill_color',[1 0 0],'fill',true,'edge',false)
plot_rf_fit(datarun_mas, glm_cellinfo.cid,'fill',true)
  
end
    
    for pair=7:12

        plot_rf_fit(datarun_mas, glm_cellinfo.pairs(1:12),'edge',true)
        plot_rf_fit(datarun_mas, glm_cellinfo.pairs(pair), 'fill_color',[1 0 0],'fill',true,'edge',false)
        plot_rf_fit(datarun_mas, glm_cellinfo.cid,'fill',true)

    end

%}

cells={6257,6256};
% XL=[0 10];

% cells={1772,1771};
% XL=[0 30];
% 
%cells={4456,4336};
% XL=[10 20];
%XL=[3 6];

celltype={'ONPar','OFFPar'};

load(['/Volumes/Analysis/nora/NSEM/BlockedSpikes/2012-08-09-3/NSEM_mapPRJ/organizedspikes_' celltype{1} '_' num2str(cells{1}) '.mat'])
S1=organizedspikes;
load(['/Volumes/Analysis/nora/NSEM/BlockedSpikes/2012-08-09-3/NSEM_mapPRJ/organizedspikes_' celltype{2} '_' num2str(cells{2}) '.mat'])
S2=organizedspikes;

Color1=[0 0.6 0.3];
Color2=[0 0.1 0.8];

hFig1=figure(1);
set(hFig1, 'Position', [100 100 1200 1000])

OP=get_cell_ids(datarun_mas,'Off Parasol');
plot_rf_fit(datarun_mas, OP(OP~=cells{2}),'fill',true,'fill_color',[0 0 0], 'edge', false, 'alpha', 0.2)
plot_rf_fit(datarun_mas, cells{1},'fill', true,'fill_color',[0 0 0], 'alpha', 0.7 , 'edge', false)
plot_rf_fit(datarun_mas, cells{2},'fill', true,'fill_color',[0 0 0], 'alpha', 0.7 , 'edge', false)
axis off
box on

%%
hFig2=figure(2);
set(hFig2, 'Position', [100 100 1400 400])

subplot(2,1,1)
hold on
for i=1:2:118
    spikes=S1.block.t_sp_withinblock{i};
    plot(spikes,i/2*ones(size(spikes)),'k.','LineWidth',0.5);
end
xlim(XL)
axis off

%%
subplot(2,1,2)
hold on
for i=1:2:118
    spikes=S2.block.t_sp_withinblock{i};
    plot(spikes,i/2*ones(size(spikes)),'k.','LineWidth',0.5);
end
xlim(XL)
axis off