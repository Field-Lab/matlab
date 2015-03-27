function [optimal_binnumber] = raster_precision(bins_per_frame, SPars, organizedspikes)
%%

bpf                 = bins_per_frame;
params.bindur       = SPars.tstim / bpf;
params.bins         = bins_per_frame *length(SPars.testframes); 
params.evalblocks   = SPars.TestBlocks;
params.trials       = length(params.evalblocks);  
params.frames       = length(SPars.testframes);
params.testdur_seconds = params.bindur * params.bins ;   
logicalspike = zeros( length(params.evalblocks) , params.bins) ;         
for i_blk = 1 : length(params.evalblocks)
	blknum = params.evalblocks(i_blk);
	sptimes = organizedspikes.block.t_sp_withinblock{blknum} - SPars.fittest_skipseconds;
	sptimes = sptimes(find(sptimes > 0 ) );
	for i_sp = 1:length(sptimes)
        spt = sptimes(i_sp);
        binnumber = ceil(spt / params.bindur );
        logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
	end
end 
clear i_blk spt sptimes  
sigmas = [.5 1 2 4 6 8 16 32];
i_sigma = 1; i_test = 1;


bins = size(logicalspike,2);

%%
bps_bysigmas = zeros(size(sigmas));
for i_sigma = 1:length(sigmas)
    sig_val        = sigmas(i_sigma);
    convolve_index = [-4*sig_val:4*sig_val];
    convolve_vec   = normpdf(convolve_index,0,sig_val) / sum(normpdf(convolve_index,0,sig_val) );
    
    
    prob_fulluop_bin   = sum(logicalspike,1) / (size(logicalspike,1));
    convolved_raw      = conv(prob_fulluop_bin, convolve_vec);
    fulluop_smooth     = convolved_raw( (4*sig_val+1):(end-4*sig_val));
    
    
    
    
    
    bps = zeros(1,5);
    for i_test = 1:5
        allreps = 1:size(logicalspike,1);
        reps     = length(allreps);
        testreps = floor(reps * .2);
        interval = floor(reps/testreps);
        testset  = (1:interval:(testreps*interval)) + i_test;
        testset  = intersect(testset, allreps);
        fitset   = setdiff(allreps, testset);
        fit_rast  = logicalspike(fitset,:);
        test_rast = logicalspike(testset,:);
        
       

        prob_uop_bin   = sum(fit_rast,1) / (size(fit_rast,1));
        prob_null_bin = (sum(prob_uop_bin) / bins) * ones(size(prob_uop_bin) );
        
        
        
        convolved_raw      = conv(prob_uop_bin, convolve_vec);
        prob_smoothuop_bin = convolved_raw( (4*sig_val+1):(end-4*sig_val));
        
        
        
        rate_smoothuop_bin  = .8* prob_smoothuop_bin + .2* fulluop_smooth; 
        rate_null_bin       = prob_null_bin;
        
        testtrials = size(test_rast,1);
        model_null      = repmat(rate_null_bin     , [testtrials,1]);
        model_smoothuop = repmat(rate_smoothuop_bin, [testtrials,1]);
        
        null_logprob            = sum(eval_rasterlogprob(test_rast, model_null, 'binary', 'conditioned'));
        smoothed_rate_logprob   = sum(eval_rasterlogprob(test_rast, model_smoothuop, 'binary', 'conditioned'));

        smoothed_rate_bits          = smoothed_rate_logprob - null_logprob;
        smoothed_rate_bits_perspike = smoothed_rate_bits / (sum(rate_null_bin ));
        bps(i_test) = smoothed_rate_bits_perspike;
        
    end
   % display(sprintf('bps for sigma of %d  equals %d', sig_val, mean(bps)));
    bps_bysigmas(i_sigma) = mean(bps);
end
sigmas = sigmas;
[blah,maxind] = max(bps_bysigmas);

optimal_binnumber = sigmas(maxind);
    
end
    
    
 


