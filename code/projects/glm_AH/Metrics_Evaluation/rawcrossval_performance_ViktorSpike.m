function [crossval_score] = rawcrossval_performance_ViktorSpike(raster_binned_sim,raster_binned_rec,Viktor_Time_Bin, bindur)
%%% PURPOSE %%%
% Compute cross-validated Fraction of Variance Explained
% Ability of Simuluated Raster to imitate (explain) recorded raster
% Only 

%%% NOTE %%%
% for really different rasters, we compute a different socre

%%% INPUTS  %%%
% bindur:= duration of bins in seconds
% Viktor_Time_Bin := 

%%% OUTPUTS %%%
% Crossval_Raw

%%% CALLS %%%
% vksp (External Viktor Spike script)

%%%%%%%%%%%%%%%
% AKHeitman 2014-12-09 


% BASIC QUANTITIES
reps        = size(raster_binned_rec,1); 
bins        = size(raster_binned_sim,2); 


% SETUP STRUCTURE
crossval_score.Viktor_Time_Bins    = Viktor_Time_Bin;
crossval_score.Viktor_Time_Secs    = bindur;
crossval_score.metric_raw          = zeros(length(Viktor_Time_Bin),1);
crossval_score.metric_name         = 'Viktor Spike';
crossval_score.normfactor_note     = 'Average raster inter_trial distance';
crossval_score.metric_note         = 'Computes Viktor Spike metric ... normalized by spike count';
crossval_score.timestamp           = datestr(clock);
crossval_score.mfile_name          = mfilename('fullpath');



%%
% LOOP THROUGH THE DIFFERENT TIME SCALES
for i_timescale = 1:length(Viktor_Time_Bin)
    viktor_cost = 1 / ( bindur * Viktor_Time_Bin(i_timescale) );
    
    % COMPUTE VIKTOR SPIKE METRIC / LOOP THROUGH  
    vksp_perspike = zeros(reps,1);
    for i_rep = 1:reps
        spt_1        = bindur * find(raster_binned_sim(i_rep,:));
        spt_2        = bindur * find(raster_binned_rec(i_rep,:));
        
        if length(spt_1) > 2*length(spt_2) || length(spt_1) < .5*length(spt_2)
            vksp_perspike(i_rep) = (length(spt_1) + length(spt_2)) /  (.5 * length(spt_1) + .5 * length(spt_2) );
            %display('Divergent Spike Trains  ');
        else
            vksp_perspike(i_rep) = spkd(spt_1, spt_2, viktor_cost) / (.5 * length(spt_1) + .5 * length(spt_2) );
        end
    end
    
    % ASSIGN
    crossval_score.metric_raw(i_timescale) = mean(vksp_perspike);
end
    
end