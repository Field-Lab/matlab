% start: 2013-09-04 AKHeitman 

% Input: Blk is just a vector of Blocks 

% Need to match exp_nm  with the mat and block matching scheme! 


function [matfiles, framenum] = blktomat_AH(exp_nm,blocks)  


fitblocks = blocks; 
blocknum  = length(fitblocks);


if strcmp(exp_nm , '2012-08-09-3') || strcmp(exp_nm, '2012-09-27-3') || strcmp(exp_nm,'2012-09-24-2') || ...
         strcmp(exp_nm, '2013-06-18-0') || strcmp(exp_nm,'2013-06-18-7') || strcmp(exp_nm, '2012-04-13-4')...
         || strcmp(exp_nm, '2013-08-19-0') || strcmp(exp_nm, '2012-09-21-3') || strcmp(exp_nm, '2013-08-19-6')

    matfiles = zeros(1, 2*blocknum);
    framenum  = zeros(1,2*blocknum*30*120);   % 2 matfiles of 30 seconds at 120 hz.
    
    for i_blk = 1:blocknum
        
        mat_ind  = [2*i_blk-1 2*i_blk];    
        stim_ind = ((i_blk-1)*(2*30*120) +1) : (i_blk*(2*30*120)) ; 

        fitblk = fitblocks(i_blk);
        startframe = (fitblk-1)*(30*120)+ 1;
        endframe   = (fitblk+1)*(30*120);

        matfiles(mat_ind) = [fitblk fitblk+1];
        framenum(stim_ind) = [startframe:endframe];
        
    end
end
    
    