% AKHeitman 2013-03-20
% Conditioned means a full matrix of rate model
% Unconditoined means a single rate model independent of trials
function [ raster_logprob_bin, logprobmat] = eval_rasterlogprob( spikemat, ratemodel, spikemat_type, ratemodel_type) 
 


trials = size(spikemat,1);
bins   = size(spikemat,2);


if strcmp(spikemat_type, 'binary')
    P_nospike = exp(-ratemodel);
    P_spike   = ratemodel.*P_nospike;
    
    spindex = find(spikemat);
    probmat = P_nospike;
    probmat(spindex) = P_spike(spindex);
    
    logprobmat         = log(probmat);
    raster_logprob_bin = mean(logprobmat,1);
end

% Could write this cleaner if I wanted to 
if ~strcmp(spikemat_type, 'binary') && strcmp(ratemodel_type,'unconditioned')
    P_nospike = exp(-ratemodel);
    P_spike   = ratemodel.*P_nospike;
    P_spike2  = (ratemodel.*P_spike ) / 2 ;
    P_spike3  = (ratemodel.*P_spike2) / 3 ;
    P_spike4  = (ratemodel.*P_spike3) / 4 ;
    P_spike5  = (ratemodel.*P_spike4) / 5 ;
    P_spike6  = (ratemodel.*P_spike5) / 6 ;
    P_spike7  = (ratemodel.*P_spike6) / 7 ;
    P_spike8  = (ratemodel.*P_spike7) / 8 ;
    P_spike9  = (ratemodel.*P_spike8) / 9 ;
    P_spike10 = (ratemodel.*P_spike9) /10 ;
    
    
	probmat = zeros(trials,bins);
    for i_trial = 1:trials
            sp_time0 = find( spikemat(i_trial ,:) == 0 );
            probmat(i_trial, sp_time0) = P_nospike(sp_time0);
    end
	for i_trial = 1:trials
            sp_time1 = find( spikemat(i_trial ,:) == 1 );
            probmat(i_trial, sp_time1) = P_spike(sp_time1);
    end
	for i_trial = 1:trials
            sp_time2 = find( spikemat(i_trial ,:) == 2 );
            probmat(i_trial, sp_time2) = P_spike2(sp_time2);
    end
	for i_trial = 1:trials
            sp_time3 = find( spikemat(i_trial ,:) == 3 );
            probmat(i_trial, sp_time3) = P_spike3(sp_time3);
    end
	for i_trial = 1:trials
            sp_time4 = find( spikemat(i_trial ,:) == 4 );
            probmat(i_trial, sp_time4) = P_spike4(sp_time4);
    end
	for i_trial = 1:trials
            sp_time5 = find( spikemat(i_trial ,:) == 5 );
            probmat(i_trial, sp_time5) = P_spike5(sp_time5);
    end
    for i_trial = 1:trials
            sp_time6 = find( spikemat(i_trial ,:) == 6 );
            probmat(i_trial, sp_time6) = P_spike6(sp_time6);
    end
    for i_trial = 1:trials
            sp_time7 = find( spikemat(i_trial ,:) == 7 );
            probmat(i_trial, sp_time7) = P_spike7(sp_time7);
    end
    for i_trial = 1:trials
            sp_time8 = find( spikemat(i_trial ,:) == 8 );
            probmat(i_trial, sp_time8) = P_spike8(sp_time8);
    end
    for i_trial = 1:trials
            sp_time9 = find( spikemat(i_trial ,:) == 9 );
            probmat(i_trial, sp_time9) = P_spike9(sp_time9);
    end
    for i_trial = 1:trials
            sp_time10 = find( spikemat(i_trial ,:) >= 10 );
            probmat(i_trial, sp_time10) = P_spike10(sp_time10);
	end
        
    
end



logprobmat  = log(probmat);
raster_logprob_bin = mean(logprobmat,1);

end

