% v18 written 2013-12-07
% Just a better variable naming than previous
% This function is largely unchanged.
% Not completely checked or noted.. but output looks correct!

%%% AK Heitman 2013-01-19
%%% seems to be working well

% Given Block Structure Input, this computes the STA
% STA computation is consistent in terms of bin placement
% .. first bin is the bin preceeding the spike!
% .. no causality issues .. just some lost time 
% .. lost time gets arbitrarily low as spike res gets better!
% Input flexibiity for different Stim Types and frame lags (frameshift)
% Flexibilty in the spike binning precision (resfactor)
% Output STA will have temporal precision of the Input Movie

%%% Perhaps have a high res output and down sampled output...


% Resfactor determines the precision of the spike binning and STA comp
%   .. number of spike bins per stim frame
%   .. integer valued  1 means spike bins equal stim frame (lowest res)
%   .. can get up to 10 (roughly millisec precision)
  
% Frameshift is a Movie effect only... no spike bins in this concept
% Accounts for possible frame delays due to hard coded Lisp Stim.
%  .ie Noise's 1st frame occurs one frame unit after initial trigger



%%% INPUTS
% blocked_stimulus: movie cell where each cell has (x,y,frames)
%        each cell corresponds to Blocks of Stimulus
% staframes:  the number of frames used in STA  (usa. 30 frames for 120 hz) 
% rawspiketimes: vector of all spike times in raw recording time units
%             has no block structure
% FitBlocks: the blocks you will use to actuall compute the STA 

% frameshift:  start at "frameshift" frames behind and go to
%          "frameshitf+staframes" behind the original spike
% resfator: max frame rate is 120 hz. but our spike recording rate is 20000 hz
%       this is a 166 fold difference
%       at resfactor of 1, we bin everything at frame rate.  ... or every  166 recording bins
%       at resfacotr of 2 we bin every half frame rate .. or every 83 collection time bins    
% resfactor should be set at 1, frameshift should be set at 1

%%% OUTPUT
% STA will be given in [x,y,t] coord.  where t is frames before the spike
% So [x,y,1] is average stim 1 frame before the spike
% [x,y,20]   is average stim 20 frames before the spike

function [STA , STAoutputnote] = STA_blockcomp_v18(blocked_stimulus,rawspiketimes,FitBlocks,blocked_frametimes,staframes,resfactor,frameshift)



n_blk = length(blocked_frametimes);
%%%
if resfactor == 1
    spikebins_starttimes = blocked_frametimes;
end

if resfactor > 1
    spikebins_starttimes = cell(n_blk,1);
    for i_blk = FitBlocks
        n_frame = length(blocked_frametimes{i_blk} ) - 1;  
        n_bin   = resfactor * n_frame;        
        spikebins_starttimes{i_blk} = linspace( blocked_frametimes{i_blk}(1) , blocked_frametimes{i_blk}(end) , (n_bin + 1) );
    end
end
block.t_sp = cell(n_blk,1); % SPIKE TIMES
block.STA = cell(n_blk,1);  % BLOCKWISE STA
block.nST = zeros(n_blk,1); 

%%
for i_blk = FitBlocks
        t_b = blocked_frametimes{i_blk}(1);
        t_e = blocked_frametimes{i_blk}(end);
        block.t_sp{i_blk} = rawspiketimes( (rawspiketimes > t_b) & (rawspiketimes <= t_e) );   
        binned_sp = nan(size(block.t_sp{i_blk})); %% just initialize vector of size t_sp
        for i_sp = 1: length(block.t_sp{i_blk})
            bin_num = max( find(spikebins_starttimes{i_blk} < block.t_sp{i_blk}(i_sp) ) );
            binned_sp(i_sp) = bin_num;
        end
     %   display(sprintf('max bin = %d' , binned_sp(end) ) )
        
        cutoff_bin = (frameshift + staframes) * resfactor ; 
        binned_sp = binned_sp ( find(binned_sp > cutoff_bin) );  %%% needs to be strictly greater
    %    remains = zeros(length(binned_sp),1);
    %    extrashift = zeros(length(binned_sp),1);
    %    standardshift = zeros(length(binned_sp),1);
    %    multiples = zeros(length(binned_sp),1);
        block.STA{i_blk} = zeros( size(blocked_stimulus{i_blk},1), size(blocked_stimulus{i_blk},2), staframes);
        clear i_sp
        for i_sp = 1:length(binned_sp);  
            %%% here use the frame shift  
          %  display(sprintf('bin is %d',  binned_sp(i_sp) ) );
            q          = binned_sp(i_sp) / resfactor ;
            spikeframe = ceil(q) - frameshift;
            remainder  = round( resfactor* ( q - floor(q) ) );
            HiRes_frontshiftweight = mod((remainder-1),resfactor)/ resfactor; 
            standardshiftweight    = (1- HiRes_frontshiftweight);
            
            standard_last  = (spikeframe-1);
            standard_first = standard_last - staframes + 1;
            
            HiRes_last      = standard_last+1;
            HiRes_first     = standard_first+1;
            
            
            SpikeMovie = standardshiftweight * blocked_stimulus{i_blk}(:,:,standard_first:standard_last)  ...
                + HiRes_frontshiftweight*blocked_stimulus{i_blk}(:,:,HiRes_first:HiRes_last);
            
            block.STA{i_blk} = block.STA{i_blk} + SpikeMovie;
        end
        block.STA{i_blk} = block.STA{i_blk} / length(binned_sp);
       % figure; imagesc(reshape(block.STA{i_blk},[size(blocked_stimulus{i_blk},1)*size(blocked_stimulus{i_blk},2) , staframes])); colorbar
        
      %  display(sprintf('done %d percent', round(100* find(FitBlocks == i_blk) / length(FitBlocks) )) ) 
end

STA = zeros(size( block.STA{FitBlocks(1)} ));
dim1 = size(STA,1) ; dim2 = size(STA,2);
for i_blk = FitBlocks
	STA = STA + block.STA{i_blk};
end
STA = STA / length(FitBlocks);

%
STAoutputnote.note1a ='STA will be given in [x,y,t] coord.  where t is frames before the spike';
STAoutputnote.note1b ='So [x,y,1] is average stim 1 frame before the spike';
STAoutputnote.note1c ='[x,y,20]   is average stim 20 frames before the spike';


STAoutputnote.note2a =  'Resfactor determines the precision of the spike binning and STA comp';
STAoutputnote.note2b =  'Resfactor is the number of spike bins per stim frame';
STAoutputnote.note2c =  'Resfactor 1 (default) means spike bins equal stim frame (lowest res)';

STAoutputnote.note3a =  'shiftfactor in terms of resolution.. where we begin STA computation';
STAoutputnote.note3b =  'Will start the STA at shiftfactor before till (shiftfactor+STAfranes)';


% Reshape to get the output into correct timing
STA2 = reshape(STA, [dim1*dim2, staframes]);
STA2 = fliplr(STA2);
STA2 = reshape(STA2, [dim1,dim2,staframes]);

STA = STA2;
          
end
        
        
        

