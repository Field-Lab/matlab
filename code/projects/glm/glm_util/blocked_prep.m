function [spikesconcat, concat_fitmovie] = blocked_prep(blockedspikes, blockedmoviecell)

% blocked spikes should be a cell of the spike times WITHIN each block
%   so each block's time should start over from 0
%
% blockedmoviecell should be cells of the blocks in the form
%   blockedmoviecell{i}.matrix = the movie in block i in (space, space,
%   time) dimensions, with a frame for each 1/120 seconds)

t_start   = 0;
tstim     = 1/120;
fitframes = size(blockedmoviecell{1}.matrix,3);
n_blocks = length(blockedspikes);

T_SP = []; blk_count = 0;
dur = tstim * fitframes;
for k = 1:n_blocks
	blk_count = blk_count + 1;
	t_sp_full = blockedspikes{k} ; % unit of time: sec, 0 for the onset of the block
	t_sp      = t_sp_full(find(t_sp_full >  t_start));
	t_sp = t_sp - t_start;
	t_spcontext = t_sp + ( blk_count -1 )*dur;
	T_SP = [T_SP ; t_spcontext];
end
spikesconcat = T_SP;

% Concatenate the fit movie (different blocks
height       = size(blockedmoviecell{1}.matrix,2);
width        = size(blockedmoviecell{1}.matrix,1);
fitframes    = size(blockedmoviecell{1}.matrix,3);
totalframes       = n_blocks * ( fitframes) ;

concat_fullfitMovie = uint8(zeros(width, height, totalframes)) ;
for i_blk = 1:n_blocks
        framenums = ( (i_blk -1)*fitframes + 1 ) :  (i_blk *fitframes);  
        concat_fullfitMovie(:,:,framenums) = blockedmoviecell{i_blk}.matrix;    
end

concat_fitmovie = concat_fullfitMovie;

end