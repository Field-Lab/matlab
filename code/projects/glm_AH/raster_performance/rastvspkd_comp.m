function [vspkd_mean vspkd_std] = rastvspkd_comp(raster_A,raster_B,timescale,bindur)
%%% PURPOSE %%%
% How far is Raster B from Raster A. Normalized by spikecount of A.
% Scores scales from 0 to 2  ( .5 is good, 1 is bad, 2 is utter failure)
%%% NOTE %%%
% for really different rasters, we compute a failure of 2
% Metric only of use within the correct range of spikes

%%% INPUTS  %%%
% raster  matrix of zeros and ones
%         rows are repitions columns are time bins
% Time scale is in seconds
% bindur  is in seconds and is the duration of raster bin

%%% OUTPUTS %%%
% Crossval_Raw

%%% CALLS %%%
% loc_spkd (External Viktor Spike script)
%          Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
%          Translated to Matlab by Daniel Reich from FORTRAN code by Jonathan Victor.

%%%%%%%%%%%%%%%
% AKHeitman 2014-12-09 

% Version 1: Try to lower input param count
% Version 0: What worked a while ago in 2014
%%

% BASIC QUANTITIES
reps        = size(raster_A,1); 
bins        = size(raster_A,2); 
cost_param  = 1/timescale;
% LOOP THROUGH THE DIFFERENT TIME SCALES

    % COMPUTE VIKTOR SPIKE METRIC / LOOP THROUGH  
vksp_perspike = zeros(reps,1);
for i_rep = 1:reps
	spt_1        = bindur * find(raster_A(i_rep,:));
	spt_2        = bindur * find(raster_B(i_rep,:));
        
    % IF SPIKE COUNT IS WAY OFF IT SHOULD BE MARKED AS A FAILURE
	if length(spt_2) > 2*length(spt_1) || length(spt_2) < .5*length(spt_1)
        vksp_perspike(i_rep) = 2;
    else
        vksp_perspike(i_rep) = loc_spkd(spt_1, spt_2, cost_param) / (length(spt_1));
	end
end
    
vspkd_mean = mean(vksp_perspike);
vspkd_std  = std(vksp_perspike);
end



function d=loc_spkd(tli,tlj,cost)
%
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