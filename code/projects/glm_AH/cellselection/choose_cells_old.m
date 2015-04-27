% AKHEITMAN 2014-11-14
% Systematic Way of Gathering Cell Lists
% Easy way of visualizing the cells which stay if we were to cut cells in a
% certain manner


clear; close all; clc
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));
selected_cells = allcells;


eval(sprintf('load %s/raster_performance/raster_scores.mat', BD.Cell_Selection));
eval(sprintf('load %s/allglmscores.mat', BD.Cell_Selection));
eval(sprintf('load %s/redo_GoodSim_Rules.mat', BD.Cell_Selection));


cut_type = 'hardrast_noglm';
if strcmp(cut_type,'softrast_noglm')
    criterion.raster.spikespersec.WN   = 5;
    criterion.raster.spikespersec.NSEM = 5;
    criterion.raster.trial_variance.WN = .25;
    criterion.raster.trial_variance.NSEM = .25;
    criterion.raster.temporal_variance.NSEM = .5;
    criterion.raster.temporal_variance.WN   = .35;
end


if strcmp(cut_type,'hardrast_noglm')
    criterion.raster.spikespersec.WN   = 10;
    criterion.raster.spikespersec.NSEM = 10;
    criterion.raster.trial_variance.WN = .2;
    criterion.raster.trial_variance.NSEM = .2;
    criterion.raster.temporal_variance.NSEM = .75;
    criterion.raster.temporal_variance.WN   = .45;
end


if strcmp(cut_type, 'norast_glm')
   criterion.glmrules = Rules{3};
end

if strcmp(cut_type, 'softrast_glm')
    criterion.glmrules = Rules{3};
    criterion.raster.spikespersec.WN   = 5;
    criterion.raster.spikespersec.NSEM = 5;
    criterion.raster.trial_variance.WN = .25;
    criterion.raster.trial_variance.NSEM = .25;
    criterion.raster.temporal_variance.NSEM = .5;
    criterion.raster.temporal_variance.WN   = .35;
end

if strcmp(cut_type,'hardrast_glm')
    criterion.glmrules = Rules{3};
    criterion.raster.spikespersec.WN   = 10;
    criterion.raster.spikespersec.NSEM = 10;
    criterion.raster.trial_variance.WN = .2;
    criterion.raster.trial_variance.NSEM = .2;
    criterion.raster.temporal_variance.NSEM = .75;
    criterion.raster.temporal_variance.WN   = .45;
end


if strcmp(cut_type, 'nocut')
    criterion = [];
end
    


i_exp = 1; i_celltype = 1;
%%


for i_exp = 1:4
    
    selected_cells{i_exp}.ONPcutcells_bycut  = [];  
	selected_cells{i_exp}.Select_ONP         = allcells{i_exp}.ONP;
	selected_cells{i_exp}.Select_ONP_index   = ones(size(allcells{i_exp}.ONP));
	selected_cells{i_exp}.Cut_ONP            = [];
    
    
    selected_cells{i_exp}.OFFPcutcells_bycut  = [];  
	selected_cells{i_exp}.Select_OFFP         = allcells{i_exp}.OFFP;
	selected_cells{i_exp}.Select_OFFP_index   = ones(size(allcells{i_exp}.OFFP));
	selected_cells{i_exp}.Cut_OFFP            = [];
  
    
    if isfield(criterion,'raster')
    for i_celltype = 1:2
        
        
        WNscores   = raster_scores{i_exp}.scores.stim_type{1}.celltype{i_celltype};
        NSEMscores = raster_scores{i_exp}.scores.stim_type{2}.celltype{i_celltype};
         
        if i_celltype == 1, cells0 = ones(size(allcells{i_exp}.ONP)); cids = allcells{i_exp}.ONP; end
        if i_celltype == 2, cells0 = ones(size(allcells{i_exp}.OFFP)); cids = allcells{i_exp}.OFFP;end
        net_cut = cells0;
        
        if isfield(criterion.raster,'spikespersec')
            WN_thresh   = criterion.raster.spikespersec.WN;
            NSEM_thresh = criterion.raster.spikespersec.NSEM;
            WN_cut = cells0; 
            WN_cut(find(WNscores.spikespersec <= WN_thresh)) = 0;
            WN_cutcells = cids(find(WNscores.spikespersec <= WN_thresh));
            NSEM_cut = cells0;
            NSEM_cut(find(NSEMscores.spikespersec <= NSEM_thresh)) = 0;
            NSEM_cutcells = cids(find(NSEMscores.spikespersec <= NSEM_thresh));
            net_cut = (NSEM_cut .* WN_cut).*net_cut;
            
            cut_cells.spikespersec.WN_cutcells  = WN_cutcells;
            cut_cells.spikespersec.NSEM_cutcells  = WN_cutcells;
       end
        
        if isfield(criterion.raster,'trial_variance')
            WN_thresh   = criterion.raster.trial_variance.WN;
            NSEM_thresh = criterion.raster.trial_variance.NSEM;
            WN_cut = cells0; 
            WN_cut(find(WNscores.trial_variance >= WN_thresh)) = 0;
            WN_cutcells = cids(find(WNscores.trial_variance >= WN_thresh));
            NSEM_cut = cells0;
            NSEM_cut(find(NSEMscores.trial_variance>= NSEM_thresh)) = 0;
            NSEM_cutcells = cids(find(NSEMscores.trial_variance >= NSEM_thresh));
            net_cut = (NSEM_cut .* WN_cut).*net_cut;
            
            cut_cells.trial_variance.WN_cutcells  = WN_cutcells;
            cut_cells.trial_variance.NSEM_cutcells  = WN_cutcells;
        end
       
        

        
   
        if i_celltype == 1 
            selected_cells{i_exp}.ONPcutcells_bycut  = cut_cells;  
            selected_cells{i_exp}.Select_ONP         = cids(find(net_cut));
            selected_cells{i_exp}.Select_ONP_index   = net_cut;
            selected_cells{i_exp}.Cut_ONP            = setdiff(cids, (cids(find(net_cut))) );
        end
        if i_celltype == 2 
            selected_cells{i_exp}.OFFPcutcells_bycut  = cut_cells;  
            selected_cells{i_exp}.Select_OFFP         = cids(find(net_cut));
            selected_cells{i_exp}.Select_OFFP_index   = net_cut;
            selected_cells{i_exp}.Cut_OFFP            = setdiff(cids, (cids(find(net_cut))) );
        end
    end
    end
end

if isfield(criterion,'glmrules');
    for i_exp = 1:4
        newcuts = criterion.glmrules.byexpnm{i_exp};
        zz = selected_cells{i_exp};
        
        zz.Select_ONP  = intersect(zz.Select_ONP, newcuts.cells);
        zz.Select_OFFP = intersect(zz.Select_OFFP, newcuts.cells);
        
        
        ONP_safeind = zeros(size(zz.ONP));
        for i_cell= 1:length(zz.Select_ONP);
            z_ind = find(zz.ONP == zz.Select_ONP(i_cell));
            ONP_safeind(z_ind ) = 1;
        end
        zz.Select_ONP_index = ONP_safeind;
        
        OFFP_safeind = zeros(size(zz.OFFP));
        for i_cell= 1:length(zz.Select_OFFP);
            z_ind = find(zz.OFFP == zz.Select_OFFP(i_cell));
            OFFP_safeind(z_ind ) = 1;
        end
        zz.Select_OFFP_index = OFFP_safeind;
        
        selected_cells{i_exp} = zz;
    end
end

        
        
        
    






%%


figure;
clf;
subplot(5,1,1);
MS = 10;
axis off
c = 0;
text(-.1, 1,sprintf('Cut criterion: %s', cut_type), 'interpreter','none');
c = c+1; text(-.1, 1-.1*c,sprintf('Xaxis: WN Normed BPS values,  Yaxis: NSEM Normed BPS,   NO SMOOTHING' ));
%c=c+1; text(-.1, 1-0.1*c,sprintf('Fit Type: %s', GLMType_Base.fitname),'interpreter','none');
c=c+1; text(-.1, 1-0.1*c,'Color are experiments, dots ONP, asterisk OFFP');
c=c+1; text(-.1, 1-0.1*c,'0 value means worse than steady firing rate,1 means unconditioned optimum');
for i_exp = 1:4
    
    if i_exp == 1; basecolor = 'r'; end
    if i_exp == 2; basecolor = 'g'; end
    if i_exp == 3; basecolor = 'b'; end
    if i_exp == 4; basecolor = 'c'; end
    
    for i_type = 1:2
        if i_type == 1
            ind = find(selected_cells{i_exp}.Select_ONP_index);
            
            
            
            
            WNvals = glmscores{i_exp}.WNBPS_ONP(ind); 
            NSEMvals = glmscores{i_exp}.NSEMBPS_ONP(ind); 
            marktype  = '.'; 
        
        end
        if i_type == 2, 
            ind = find(selected_cells{i_exp}.Select_OFFP_index);
            WNvals = glmscores{i_exp}.WNBPS_OFFP(ind);  
            NSEMvals = glmscores{i_exp}.NSEMBPS_OFFP(ind);
            marktype = '*'; 
        end
        
        subplot(5,2, (i_exp*2 + i_type))
        hold on;
        WNvals( WNvals <=0) = 0;
        NSEMvals( NSEMvals <=0) = 0;
        
        if ~isempty(WNvals), max_x = max(1, max(WNvals)); end
        if ~isempty(NSEMvals),max_y = max(1, max(NSEMvals)); end
        max_val = max(max_x,max_y);
        
        set(gca,'xlim',[0 max_x]); set(gca,'ylim',[0 max_y]);
        
        
        plotstring = sprintf('%s%s',basecolor,marktype);
        plot(WNvals, NSEMvals, plotstring,'markersize',MS);
        plot(linspace(0,1,100), linspace(0,1,100), 'k')
    end
end
orient tall
eval(sprintf('print -dpdf %s/WN_vs_NSEM_%s', BD.Cell_Selection, cut_type))

            