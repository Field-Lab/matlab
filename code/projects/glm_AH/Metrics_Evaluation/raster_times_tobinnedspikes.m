
function binned_sp = raster_times_tobinnedspikes( sptimecell , bindur , timeinterval);

ta = timeinterval(1);
tb = timeinterval(2);
tdur = tb - ta;
bins = ceil( tdur / bindur); 

trials = length(sptimecell);


testbins = bpf *length(SPars.testframes);   
trials = length(evalblocks);       
        logicalspike = zeros( length(evalblocks) , testbins) ; 
        
        for i_blk = 1 : length(evalblocks)
            blknum = evalblocks(i_blk);
            sptimes = organizedspikes.block.t_sp_withinblock{blknum} - SPars.testseconds(1);
            sptimes = sptimes(find(sptimes > 0 ) );
            for i_sp = 1:length(sptimes)
                spt = sptimes(i_sp);
                binnumber = ceil(spt / bindur );
                logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
            end
        end



end