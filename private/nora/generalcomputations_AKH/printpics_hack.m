clear; close all

%%
%stimname = 'eye-long-v2.rawMovie';
%NSEMmoviedir = sprintf('/netapp/snle/lab/Experiments/Array/Analysis/akheitman/Decomposition/old_NSEM_directsaving_%s', stimname);


%stimname = 'eye-120-3_0-3600.rawMovie';
%NSEMmoviedir = sprintf('/netapp/snle/lab/Experiments/Array/Analysis/akheitman/Decomposition/NSEM_directsaving_%s', stimname);


stimname = 'eye-120-3_0-3600.rawMovie';
NSEMmoviedir = sprintf('/netapp/snle/lab/Experiments/Array/Analysis/akheitman/Decomposition/NSEM_Movie');

picdir = sprintf('%s/pics3',NSEMmoviedir);
if ~exist(picdir, 'dir'), mkdir(picdir); end
%savedir = sprintf('/netapp/snle/lab/Experiments/Array/Analysis/akheitman/Decomposition/NSEM_ConcatinTime_SplitbyColumn_%s', stimname);
%if ~exist(savedir), mkdir(savedir); end
%framesperblock = 3600;
figure;
for i_block = 1:2
    
    eval(sprintf('load %s/eyemov%d.mat X Xr X_rescaled', NSEMmoviedir, i_block))
    
    for i_pic = 1:30
        frame = 120* i_pic;
        clf;
        subplot(3,1,1)
        imagesc((squeeze(X(frame,:,:)))'); colormap gray
        subplot(3,1,2)
        imagesc((squeeze(Xr(frame,:,:)))'); colormap gray
        subplot(3,1,3)
        imagesc((squeeze(X_rescaled(frame,:,:)))'); colormap gray
        
        imagenum = (i_block-1)*30 + i_pic;
        orient tall
        eval(sprintf('print %s/image%d.pdf -dpdf', picdir, imagenum))
    end
    
        
    

end
        
    