% AKHeitman 2014-04-14

%%% PURPOSE %%%
% Concatenate the fit movie into a single long matrix



% just to keep my fitting consistent
function concat_fitmovie = concat_fitmovie_fromblockedcell(blockedmoviecell , FitPars)

% FitPars needs
%   .width
%   .height
%   .FitBlocks
%   .novelblocks
%   .fitframes

height       = FitPars.height;
width        = FitPars.width;
fitblocks    = FitPars.FitBlocks;
fitframes    = FitPars.fitframes;
novelblocks  = FitPars.NovelBlocks;

fitframesperblock = length(fitframes) ;
totalframes       = length(fitblocks) * ( fitframesperblock) ;
concat_fullfitMovie = uint8(zeros(width, height, totalframes)) ;
for i_blk = fitblocks
        blkind = find(fitblocks == i_blk);
        framenums = ( (blkind -1)*fitframesperblock + 1 ) :  (blkind *fitframesperblock);  
        n_blkind = find(novelblocks == i_blk);
        concat_fullfitMovie(:,:,framenums) = blockedmoviecell{n_blkind}.matrix (:,:, fitframes);    
end

concat_fitmovie = concat_fullfitMovie;

end