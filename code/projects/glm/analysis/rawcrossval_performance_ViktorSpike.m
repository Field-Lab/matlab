function [crossval_score] = rawcrossval_performance_ViktorSpike(raster_binned_sim,raster_binned_rec,Viktor_Time_Secs, bindur)
%%% PURPOSE %%%
% Compute cross-validated Fraction of Variance Explained
% Ability of Simuluated Raster to imitate (explain) recorded raster
% Only 

%%% NOTE %%%
% for really different rasters, we compute a different socre

%%% INPUTS  %%%
% bindur:= duration of bins in seconds
% Viktor_Time_Bin

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
crossval_score.Viktor_Time_Secs    = Viktor_Time_Secs;
crossval_score.metric_raw          = zeros(length(Viktor_Time_Secs),1);
crossval_score.metric_name         = 'Viktor Spike';
crossval_score.normfactor_note     = 'Average raster inter_trial distance';
crossval_score.metric_note         = 'Computes Viktor Spike metric ... normalized by spike count';
crossval_score.timestamp           = datestr(clock);
crossval_score.mfile_name          = mfilename('fullpath');



%%
% LOOP THROUGH THE DIFFERENT TIME SCALES
for i_timescale = 1:length(Viktor_Time_Secs)
    viktor_cost = 1 / ( Viktor_Time_Secs );    
    % COMPUTE VIKTOR SPIKE METRIC / LOOP THROUGH  
    vksp_perspike = zeros(reps,1);
    for i_rep = 1:reps
        spt_1        = bindur * find(raster_binned_sim(i_rep,:));
        spt_2        = bindur * find(raster_binned_rec(i_rep,:));
        
        % Special Caveat for when the simulation blows up
        % Only divide by spikes of recorded .. how far is sim from recorded
        if length(spt_1) > 2*length(spt_2)
            vksp_perspike(i_rep) = (length(spt_1)) /  (length(spt_2));
            %display('Divergent Spike Trains  ');
        else
            vksp_perspike(i_rep) = spkd(spt_1, spt_2, viktor_cost) / (length(spt_2) );
        end
    end    
    % ASSIGN
    crossval_score.metric_raw(i_timescale) = mean(vksp_perspike);
    crossval_score.mean_vspkd(i_timescale) = mean(vksp_perspike);
    crossval_score.std_vspkd(i_timescale)  = std(vksp_perspike);
    crossval_score.compared_spiketrains_count = reps;
end
    
end


function d=spkd(tli,tlj,cost)
% Copied from http://www-users.med.cornell.edu/~jdvicto/spkdm.html
% 2014-11-17


% d=spkd(tli,tlj,cost) calculates the "spike time" distance
% (Victor & Purpura 1996) for a single cost
%
% tli: vector of spike times for first spike train
% tlj: vector of spike times for second spike train
% cost: cost per unit time to move a spike
%
%  Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
%  Translated to Matlab by Daniel Reich from FORTRAN code by Jonathan Victor.
%
nspi=length(tli);
nspj=length(tlj);

if cost==0
   d=abs(nspi-nspj);
   return
elseif cost==Inf
   d=nspi+nspj;
   return
end

scr=zeros(nspi+1,nspj+1);
%
%     INITIALIZE MARGINS WITH COST OF ADDING A SPIKE
%
scr(:,1)=(0:nspi)';
scr(1,:)=(0:nspj);
if nspi & nspj
   for i=2:nspi+1
      for j=2:nspj+1
         scr(i,j)=min([scr(i-1,j)+1 scr(i,j-1)+1 scr(i-1,j-1)+cost*abs(tli(i-1)-tlj(j-1))]);
      end
   end
end
d=scr(nspi+1,nspj+1);

end

