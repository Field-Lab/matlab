function [raster_binned] = blockedspiketimes_to_binnedraster(blocked_spiketimes,raster_params)

%%% PURPOSE %%%
% Take raw spike times in blocks (time relative to block start)
% Spit out a matrix raster binned according to raster_params


%%% INPUTS  %%%
% Blocked_SpikeTimes
%  usa. taken from organizedspikes.blocks (seperate function written AKHEITMAN)
%  Spike times in each recording block
%  Not absolute time, but time relative beginning of each block
% Raster_Params
%  .bins := number of bins needed
%  .bindur := duration of bins in seconds
%  .evalblocks := which blocks to draw spike times froms
%  .fittest_skipseconds := how many seconds to ignore up front (adaptation
%  to gray screen or changes in stimulus structure)

%%% OUTPUTS %%%
% Raster_Binned
%  [trials,bins] matrix of zeros and ones
%  rows are different trials (blocks of spikes)
%  bins is just time bins
%  zeros = no spike , one = spike

%%% CALLS %%%
% none


% AKHEITMAN 2014-12-04

raster_binned = zeros( length(raster_params.evalblocks) , raster_params.bins) ;         
for i_blk = 1 : length(raster_params.evalblocks)
    blknum = raster_params.evalblocks(i_blk);
    sptimes = blocked_spiketimes{blknum} - raster_params.fittest_skipseconds;
    sptimes = sptimes(find(sptimes > 0 ) );
    for i_sp = 1:length(sptimes)
        spt = sptimes(i_sp);
        binnumber = ceil(spt / raster_params.bindur );
        if binnumber <= raster_params.bins
            raster_binned( i_blk, binnumber )  =  raster_binned( i_blk,binnumber ) + 1;
        end
    end
end

end
                