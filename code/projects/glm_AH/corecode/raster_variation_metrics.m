function [Raster_Variability] = raster_variation_metrics( RasterCell, raster_time, start_delay, end_early, maxpairs)
%%% AK Heitman
%%% checking on 2012-12-06  Good Complete and Easy to read!

% Computes pairwise spike metric distances of a given raster
% prints out summart stats of on avg how varaible the rasters are.. running
% throuh a 5 different spke metric params  (millisec, centisec, decisec,
% and in between )

% Calls: spkd  by Victor to compute individual comparisons

% Does the minimum of either all pairs or maxpairs
% RasterCell is a cell with trials numbers of entry.. each
% RasterCell{i}.spikes  i ranges from 1 to trial number  where .spikes is
% just a vector of spike times

% start_delay : how much of the front raster we want to discount in seconds
% end_early : how much of the end of raster we want to discount in seconds
% maxpairs: the number of pairs blows up, so this puts a cap on how much we
% compute.. shouldn't tak e too many to converge to an average distance



% Spits out RasterVariability   which shows all pair data in
% .individualpairs
% more important is .summary_stats with the avg distances
%.start_Delat and .end_early are the input parameters



%%%%%%%% LOAD AND FORMAT SPIKES FOR METRIC CALCULATION 
trials = max (size(RasterCell) );  
retsp_formetric = cell (trials , 1);
for i_trial = 1 : trials
    ret_spike        = RasterCell{i_trial}.spikes; 
	ret_ind1         = find(ret_spike > start_delay ) ;
	ret_ind2         = find(ret_spike < raster_time - end_early);
	ret_index        = intersect(ret_ind1 , ret_ind2);
	ret_final        = ret_spike (ret_index );        
    retsp_formetric{i_trial, 1}.metricspikes      = ret_final;
	retsp_formetric{i_trial, 1}.trial             = i_trial; 
    
	retsp_formetric{i_trial, 1}.originalspikesecs = ret_spike;
end


%%%%%%%% CHOOSE PAIRS FOR METRIC COMPUTATION
allpairs       = nchoosek( 1:trials , 2);
totalpairs     = size(allpairs,1);
if ~exist('maxpairs' ,'var')
    maxpairs  = totalpairs + 1;
end
talliedpairs   = min(maxpairs,totalpairs); 
if talliedpairs < totalpairs
    chosenones     = ceil ( totalpairs * rand (talliedpairs , 1) ) ;
else
    chosenones     = 1:totalpairs;
end
chosenpairs    = allpairs(chosenones,:);


%%%%%%%% COMPUTE METRICS BETWEEN PAIRS .. 5 DIFFERENT TIME SCALES
Pairwise.pairs = chosenpairs;
%Pairwise.metricdist_millisec     = zeros ( talliedpairs, 1);
Pairwise.metricdist_halfcentisec = zeros ( talliedpairs, 1);
Pairwise.metricdist_centisec     = zeros ( talliedpairs, 1);
Pairwise.metricdist_halfdecisec  = zeros ( talliedpairs, 1);
%Pairwise.metricdist_decisec      = zeros ( talliedpairs, 1);

%%%%%%%% COMPUTE METRICS AMONG PAIRS
pctdonepre = 0;
for i_pair = 1 : talliedpairs
	train1 = retsp_formetric{chosenpairs(i_pair , 1) , 1}.metricspikes;
	train2 = retsp_formetric{chosenpairs(i_pair , 2) , 1}.metricspikes;
%	Pairwise.metricdist_millisec(i_pair)    = spkd ( train1 , train2 , 1000 );
	Pairwise.metricdist_halfcentisec(i_pair)= spkd ( train1 , train2 , 200  );
	Pairwise.metricdist_centisec(i_pair)    = spkd ( train1 , train2 , 100  );
	Pairwise.metricdist_halfdecisec(i_pair) = spkd ( train1 , train2 , 20   );
%	Pairwise.metricdist_decisec(i_pair)     = spkd ( train1 , train2 , 10   );  
	pctdone =100* i_pair / talliedpairs;
	if round(pctdone/10) > round(pctdonepre/10);
        display(sprintf('PercentRasterMetricVariabilityDone_%d',round(pctdone)));
	end
	pctdonepre = pctdone;
end
clear pctdonepre pctdone

%%%%%%%% SUMMARY STATISTICS
%ret_metricvar.avgmillisec        = mean(Pairwise.metricdist_millisec);
%ret_metricvar.stdmillisec        =  std(Pairwise.metricdist_millisec);
ret_metricvar.avghalfcentisec    = mean(Pairwise.metricdist_halfcentisec);
ret_metricvar.stdhalfcentisec    =  std(Pairwise.metricdist_halfcentisec);
ret_metricvar.avgcentisec        = mean(Pairwise.metricdist_centisec);
ret_metricvar.stdcentisec        =  std(Pairwise.metricdist_centisec);
ret_metricvar.avghalfdecisec     = mean(Pairwise.metricdist_halfdecisec);
ret_metricvar.stdhalfdecisec     =  std(Pairwise.metricdist_halfdecisec);
%ret_metricvar.avgdecisec         = mean(Pairwise.metricdist_decisec);
%ret_metricvar.stddecisec         =  std(Pairwise.metricdist_decisec); 

%%%%%%%% ASSIGN THE OUTPUT VARIABLES
Raster_Variability.summary_stats        = ret_metricvar;
Raster_Variability.individualpairs      = Pairwise;
Raster_Variability.start_delay          = start_delay;
Raster_Variability.end_early            = end_early;
Raster_Variability.duration             = raster_time;
Raster_Variability.trials               = trials;
Raster_Variability.pairstested          = talliedpairs;