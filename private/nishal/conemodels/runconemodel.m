% AKHeitman 2013-01-26
% Runs cone model over a spatial temporal movie (space, space, time)
% One pixel at a time.
% Assumes an input of Rstar_frame .. output is in pAmp;

% THESE MODELS SEEM TO RUN MUCH FASTER ON BERTHA!! 
%%% Test Call %%%
%{



%}

function pAmpmovie = runconemodel(model_name, movie_rstar_frame, framedur, binsperframe, checkfunction)



%%%
%%
if strcmp(model_name, 'model1')
    path('./model1_Juan_20140122', path)
    load('./model1_Juan_20140122/EyeMovements_092413Fc12.mat')
    load('./model1_Juan_20140122/DimFlash_092413Fc12.mat')
end
pAmpmovie = zeros(size(movie_rstar_frame));

columns =size(movie_rstar_frame,1);
rows   = size(movie_rstar_frame,2);
frames = size(movie_rstar_frame,3);
time_frame  = framedur : framedur : (frames * framedur);
time_bin  = (framedur/binsperframe): (framedur/binsperframe): (frames * framedur);

%i_row = 1; i_col = 1;
%%
for i_row = 1:rows
    display(sprintf('Working on Row %d of %d', i_row, rows));
    for i_col = 1:columns
%%
        pixmovie_perframe = squeeze(movie_rstar_frame(i_col,i_row,:))';
        pixmovie_perbin = (1/binsperframe)*reshape( repmat(pixmovie_perframe, binsperframe,1) , binsperframe*length(pixmovie_perframe),1) ;
        pixmovie_perbin = pixmovie_perbin';
        
        if i_row == 1 && i_col == 1
            if exist('checkfunction') && checkfunction.check
                clf; subplot(2,1,1); plot(time_frame,pixmovie_perframe); hold on; title('Check upsampling worked'); ylabel('R*perframe'); hold off
                subplot(2,1,2); plot(time_bin, pixmovie_perbin); hold on; ylabel('R*perbin'); xlabel('seconds'); 
                filename= sprintf('%s_CheckingUpsamlingandRstarcounts.pdf',checkfunction.plotblockname);
                orient tall
                eval(sprintf('print -dpdf %s/%s', checkfunction.plotdir, filename));
            end
        end
        %%
        %%%%% RUN THE MODEL %%%%%%%%%
        switch model_name
            case {'model1'}
                cmodel_output=-riekemodelstim_slowadapt([DimFlash.coeffs(5) EyeMovements.slowcacoeff],time_bin,pixmovie_perbin, DimFlash.coeffs(1:7));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i_row == 1 && i_col == 1
            if exist('checkfunction') && checkfunction.check
                stimnormed = (pixmovie_perbin - min(pixmovie_perbin)) / max(pixmovie_perbin - min(pixmovie_perbin));
                cmodel_outputnormed  =  (cmodel_output - min(cmodel_output)) / max(cmodel_output - min(cmodel_output));
                origstim_normed = (pixmovie_perframe - min(pixmovie_perframe)) / max(pixmovie_perframe - min(pixmovie_perframe));
                cmodel_outputnormed_ds = interp1(time_bin,cmodel_outputnormed,time_frame); 
                clf;
                subplot(4,1,1); plot( time_bin ,  pixmovie_perbin, 'b'); ylabel('R*perbin')
                subplot(4,1,2); plot( time_bin , cmodel_output, 'r'); ylabel('pA')
                subplot(4,1,3); hold on; plot(time_bin ,stimnormed,'b'); plot(time_bin,cmodel_outputnormed, 'r'); xlabel('Original output in bins');
                subplot(4,1,4); hold on; plot(time_frame ,origstim_normed,'b');ylabel('normed'); plot(time_frame,cmodel_outputnormed_ds, 'r'); xlabel('Downsampled to frames');
                
                filename = sprintf('%s_cmodeloutput_orgistim',checkfunction.plotblockname);
                orient tall
                eval(sprintf('print -dpdf %s/%s', checkfunction.plotdir, filename));
                
            end
        end
        cmodel_frame = interp1(time_bin,cmodel_output,time_frame);        
        pAmpmovie(i_col,i_row,:) = cmodel_frame;
    end
end


end
