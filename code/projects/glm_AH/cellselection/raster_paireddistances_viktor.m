function [paired_distance] = raster_paireddistances_viktor(raster_binned,metric_params)
%%% PURPOSE %%%
% Compute inter-repition distance in raster trials
% Compute distance Viktor Spikes, with various time parameters
% Average distance over a number of pairs

%%% INPUTS  %%%
% Raster_Binned
%  [trials,bins] matrix of zeros and ones
%  rows are different trials (blocks of spikes)
%  bins is just time bins
%  zeros = no spike , one = spike
% Metric_params
%  .ViktorTimeParams_Bins: bins scale of the Viktor cost parameter
%                cost parameter computed as (1/x), it costs one to move x
%  .bindur:  duration of each bin in seconds (final input into vksp)
%  .pairnumbers: number of random pairs used to compute avg paired distance

%%% OUTPUTS %%%
% Paired_distance


%%% DETAILS %%%
% Normalize Voktor SPike disance by the avg number of spikes

%%% MAJOR CALLS %%%
% spkd : The core viktor spike metric, external code

% AKHeitman 2014-12-05
% AKHeitman 2014-12-08 clean up / comment

%%
% UNPACK PARAMETERS
VikTime     = metric_params.ViktorTimeParams_Bins;
bindur      = metric_params.bindur;
pairnumbers = metric_params.pairnumbers;

% INTIAILIZE THE SOLUTION STRUCTURE 
clear paired_distance
paired_distance.Viktordist          = zeros(length(VikTime),1);
paired_distance.Viktordist_perspike = zeros(length(VikTime),1);
paired_distance.ViktorParams_note    = sprintf('Cost computed with cost (1/(x)) where x is time of smoothing var with %d pairs' , pairnumbers);
paired_distance.ViktorParams_Bins   = VikTime;
paired_distance.ViktorParams_Time   = bindur*VikTime;
paired_distance.timestamp           = datestr(clock);
paired_distance.mfile_name          = mfilename('fullpath');

% LOCAL PARAMS
pars.bins           = size(raster_binned,2);
pars.reps           = size(raster_binned,1);
pars.seconds        = bindur *size(raster_binned,2);
pars.spikes_pertrial= sum(raster_binned(:)) / pars.reps;
pars.spikes_persec  = pars.spikes_pertrial / ( bindur * size(raster_binned,2) );

% LOOP THROUGH TIMESCALES
for i_timescale = 1:length(VikTime);
    
    % COMPUTE VIKTOR SPIKE METRIC / LOOP THROUGH
    vksp    = zeros(pairnumbers,1);    
    vksp_perspike = zeros(pairnumbers,1);
    viktor_cost = 1 / ( bindur*VikTime(i_timescale) );
    for i_pair = 1:pairnumbers 
        pair_vals = randsample(pars.reps,2);
        ind1 = pair_vals(1);
        ind2 = pair_vals(2);
        spt_1        = bindur * find(raster_binned(ind1,:));
        spt_2        = bindur * find(raster_binned(ind2,:));
        vksp(i_pair) = spkd(spt_1, spt_2, viktor_cost);
        vksp_perspike(i_pair) = spkd(spt_1, spt_2, viktor_cost) / (.5 * length(spt_1) + .5 * length(spt_2) );
    end
    paired_distance.Viktordist(i_timescale)                         = mean(vksp) ; 
    paired_distance.Viktordist_perspike(i_timescale)                = mean(vksp_perspike) ; 
    
end

end