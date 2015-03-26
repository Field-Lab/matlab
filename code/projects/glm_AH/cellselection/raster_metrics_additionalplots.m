% Designed to further explore the raster metrics with respect to whatever


clear; close all; clc
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));  % allcells structire
expanded_metrics = true;
rastmetdir = sprintf('%s/raster_performance', BD.Cell_Selection);
if ~exist(rastmetdir, 'dir'), mkdir(rastmetdir); end

raster_scores = allcells;

DirPars.rast_dir = BD.BlockedSpikes;

params.bindur     = .00083275;
params.bins_per_frame = 10;
metricparams.smoothbins    = 10;
metricparams.bindur = params.bindur;
if exist('expanded_metrics', 'var') && expanded_metrics
    metricparams.slowsmoothbins    = 25; % roughly 20 msecs
    metricparams.slowestsmoothbins = 75; %roughly 50 msecs 
    metricparams.pair_numbers      = 40;
end
    


eval(sprintf('load %s/expanded_raster_scoresA.mat', rastmetdir));
eval(sprintf('load %s/expanded_raster_scoresB.mat', rastmetdir));
eval(sprintf('load %s/expanded_raster_scoresC.mat', rastmetdir));
eval(sprintf('load %s/expanded_raster_scoresD.mat', rastmetdir));


eval(sprintf('load %s/rastA_byeye.mat', rastmetdir));
eval(sprintf('load %s/rastB_byeye.mat', rastmetdir));
eval(sprintf('load %s/rastC_byeye.mat', rastmetdir));
eval(sprintf('load %s/rastD_byeye.mat', rastmetdir));

close all


%%

i_metric = 6 ; i_time = 1; i_exp = 1; i_type = 1;
for i_metric = 1:6
    for i_time = 1:3
        
        if i_metric == 1, title_base = 'Temporal Variance Normed';  pdf_title = 'TemporalVariance';     end
        if i_metric == 2, title_base = 'Bits Per Spike Mean';       pdf_title = 'BPS_Mean';             end
        if i_metric == 3, title_base = 'Bits Per Spike Max';        pdf_title = 'BPS_Max';              end
        if i_metric == 4, title_base = 'Bits Per Spike Std';        pdf_title = 'BPS_Std';              end
        if i_metric == 5, title_base = 'Fractional Error Capped at 30';          pdf_title = 'pair_distance_vector'; end
        if i_metric == 6, title_base = 'Viktor Distance Per Spike'; pdf_title = 'pair_distance_Viktor'; end
        
        if i_time ==1 , bins = metricparams.smoothbins; end
        if i_time ==2 , bins = metricparams.slowsmoothbins; end
        if i_time ==3 , bins = metricparams.slowestsmoothbins; end
        


        MS = 12;
        %%
        for i_exp = 1:4
            clf;
            if i_exp == 1, expname = 'expA'; rast_byeye = rastA_byeye; fullscores0 = expanded_raster_scoresA; end
            if i_exp == 2, expname = 'expB'; rast_byeye = rastB_byeye; fullscores0 = expanded_raster_scoresB; end
            if i_exp == 3, expname = 'expC'; rast_byeye = rastC_byeye; fullscores0 = expanded_raster_scoresC; end
            if i_exp == 4, expname = 'expD'; rast_byeye = rastD_byeye; fullscores0 = expanded_raster_scoresD; end
            
            fullscores = fullscores0.expand{i_time};
            ONP  = allcells{i_exp}.ONP;
            OFFP = allcells{i_exp}.OFFP;
            
            celllistONOFF = rast_byeye(:,1);
            for i_type = 1:2
              
                goodWNcells     = celllistONOFF(intersect( find(rast_byeye(:,(i_type+1))) , find(rast_byeye(:,4))  ));
                goodNSEMcells   = celllistONOFF(intersect( find(rast_byeye(:,(i_type+1))) , find(rast_byeye(:,6))  ));
                badWNcells      = celllistONOFF(intersect( find(rast_byeye(:,(i_type+1))) , find(rast_byeye(:,5))  ));
                badNSEMcells    = celllistONOFF(intersect( find(rast_byeye(:,(i_type+1))) , find(rast_byeye(:,7))  ));
                
                if i_type == 1, cells = ONP;  marktype = '.'; ctype = 'ONP'; end
                if i_type == 2, cells = OFFP; marktype = '*'; ctype = 'OFFP';end
                
                wn_struct   = fullscores.stim_type{1}.celltype{i_type};
                nsem_struct = fullscores.stim_type{2}.celltype{i_type};
                if i_metric == 1, wn_vec = sqrt(wn_struct.temporal_variance_normed); end
                if i_metric == 2, wn_vec = wn_struct.bps_mean;    end
                if i_metric == 3, wn_vec = wn_struct.bps_max; end
                if i_metric == 4, wn_vec = wn_struct.bps_std;    end
                if i_metric == 5, wn_vec = wn_struct.fractional_error;    end
                if i_metric == 6, wn_vec = wn_struct.viktordist_perspike;    end

                if i_metric == 1, nsem_vec = sqrt(nsem_struct.temporal_variance_normed); end
                if i_metric == 2, nsem_vec = nsem_struct.bps_mean;    end
                if i_metric == 3, nsem_vec = nsem_struct.bps_max; end
                if i_metric == 4, nsem_vec = nsem_struct.bps_std;    end
                if i_metric == 5, nsem_vec = nsem_struct.fractional_error;    end
                if i_metric == 6, nsem_vec = nsem_struct.viktordist_perspike;    end
                if i_metric == 5, wn_vec(wn_vec>30) = 30; nsem_vec(nsem_vec>30)=30; end
                
                for i_col = 1:2
                    if i_col == 1, goodcells = goodWNcells;     badcells = badWNcells;      cut_type = 'WN';     end
                    if i_col == 2, goodcells = goodNSEMcells;   badcells = badNSEMcells;    cut_type = 'NSEM';   end 
                    
                    goodcells_index = [];
                    for i_cell = 1:length(goodcells); 
                        goodcells_index= [goodcells_index; find(cells == goodcells(i_cell))]; 
                    end
                    badcells_index = [];
                    for i_cell = 1:length(badcells); 
                        badcells_index= [badcells_index; find(cells == badcells(i_cell))]; 
                    end
                    
                    subplot(2,2,2*(i_type-1) + i_col)
                    hold on; set(gca,'fontsize', 10)
                    timeadded = sprintf('Time %dmsecs', round(1000*params.bindur*bins) );
                    pdftime = sprintf('timescale_%dBins', bins);
                    title(sprintf('%s %s, Cut by: %s,Cells: %s', title_base, timeadded,cut_type,ctype));
                    xlabel('WN values');
                    ylabel('NSEM values');
                    
                    standardstring = sprintf('k%s',marktype);
                    goodstring = sprintf('g%s',marktype);
                    badstring = sprintf('r%s',marktype);

                    plot(wn_vec, nsem_vec, standardstring,'markersize',MS);  
                    plot(wn_vec(goodcells_index), nsem_vec(goodcells_index), goodstring,'markersize',MS); 
                    plot(wn_vec(badcells_index), nsem_vec(badcells_index), badstring,'markersize',MS); 
                    
                    max_x = 0;
                    max_y = 0;
                    max_x = max(max_x, max(wn_vec));
                    max_y = max(max_y, max(nsem_vec)); 
                    max_val = max(max_x, max_y);
                    
                    set(gca,'xlim',[0 max_val]); set(gca,'ylim',[0 max_val]);
                    plot(linspace(0,max_val,100), linspace(0,max_val,100),'k');
                end
            end
            
            orient landscape
            eval(sprintf('print -dpdf %s/%s_%s_%s.pdf',rastmetdir, pdf_title,pdftime, expname)); 
        end
        %%
        
        
        
        
    end
end