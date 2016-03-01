% AKHeitman 2014-04-27
% GLMPars = GLMParams  or GLMPars = GLMParams(GLMType.specialchange_name)
% Sensitive to naming changes in GLMParams.
% Only saves stuff (and calls directories) if we are in a troubleshooting mode
% Heavily GLMType dependent computations will be carried out here
% Outsourcable computations will be made into their own functions
%troubleshoot optional
% need troublshoot.doit (true or false), 
%troubleshoot.plotdir,
% troubleshoot.name
function [results] = prep_spikefiltersGP(GLMType, spikes, fitmovie, center_coord,STA,troubleshoot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up GLMParams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GLMPars = GLMParams;
if GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end
if GLMType.debug, GLMPars.tolfun = 2;  end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cut the Movie down to ROI
% Normalize Movie and set nullpoint
% output of this section is stim
% stim in [xy, time] coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROI_length      = GLMPars.stimfilter.ROI_length;
stimsize.width  = size(fitmovie,1);
stimsize.height = size(fitmovie,2);
stimsize.frames = size(fitmovie,3);
ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
stim            = fitmovie(ROIcoord.xvals, ROIcoord.yvals, :);


fitmoviestats.mean     =  double(mean(fitmovie(:)));
fitmoviestats.minval   =  double(min(fitmovie(:)));
fitmoviestats.maxval   =  double(max(fitmovie(:)));
fitmoviestats.span     =  fitmoviestats.maxval - fitmoviestats.minval;
fitmoviestats.normmean =  fitmoviestats.mean   / fitmoviestats.span ;

stim   = double(stim);
stim   = stim - fitmoviestats.minval;
stim   = stim / fitmoviestats.span;

if strcmp(GLMType.nullpoint, 'mean')
    stim = stim - fitmoviestats.normmean;
else
    error('you need to fill in how to account for stimulus with a different nullpoint')
end

if exist('troubleshoot','var') && troubleshoot.doit
    clf
    subplot(3,2,[3 5]); hist(double(fitmovie(:)),20); set(gca, 'fontsize', 10); title('histogram of full raw stim')
    subplot(3,2,[4 6]); hist(stim(:) ,20); set(gca, 'fontsize', 10); title('histogram of normed stim over ROI') 
    subplot(3,2,[1 2]); set(gca, 'fontsize', 10); axis off
    c = 0;
    c=c+1; text(-.1, 1-0.1*c,sprintf('Trouble shooting: %s',troubleshoot.name));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Specifically: Stimulus Normalizing Component of "glm-nospace"' ));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Plot Date %s',datestr(clock)));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')) );
    
    orient landscape
    eval(sprintf('print -dpdf %s/%s_glmnospace_stimnorm.pdf', troubleshoot.plotdir, troubleshoot.name));
end 
stim = reshape(stim, [ROI_length^2 , stimsize.frames]);

clear stimsize ROIcoord ROI_length fitmoviestats 

%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create spfilter from the STA 
% output: spfilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROI_length      = GLMPars.stimfilter.ROI_length;
stimsize.width  = size(fitmovie,1);
stimsize.height = size(fitmovie,2); 
ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
spfilter = spatialfilterfromSTA(STA,ROIcoord.xvals,ROIcoord.yvals);
if exist('troubleshoot','var') && troubleshoot.doit
    clf
    subplot(3,2,[1 2]);   set(gca, 'fontsize', 10); axis off
    c = 0;
    c=c+1; text(-.1, 1-0.1*c,sprintf('Trouble shooting: %s',troubleshoot.name));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Specifically: spfilter from WN-STA "glm-nospace/spatialfilterfromSTA"' ));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Plot Date %s',datestr(clock)));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')) );
    
    
    subplot(3,2,[3 5]);  set(gca, 'fontsize', 10); imagesc(reshape(spfilter,[ROI_length, ROI_length])); colorbar
    xlabel('pixels'); ylabel('pixels');
    title('Spatial Filter first rank of the STA');
    
    subplot(3,2,[4 6]); set(gca, 'fontsize', 10);   imagesc( (squeeze(mean(STA,1))' )); colorbar
    ylabel('frames');  xlabel('pixel rows')
    title('Raw STA, columns collapsed to 1-spatial dimension');
    
    
    orient landscape
    eval(sprintf('print -dpdf %s/%s_glmnospace_spfilterfromSTA.pdf', troubleshoot.plotdir, troubleshoot.name));
end 
spfilter  = spfilter';  % dimension [1,ROI_length]

clear stimsize ROIcoord ROI_length fitmoviestats 

%%
%%%%%%%%%%%%%%%%%%%%%%%
% Creat the final stim dependent input (X) to the optimization algorithm
% X_bin  ([rank, bins])
%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(GLMType.k_filtermode, 'fixedSP_rk1_linear')
    if strcmp(GLMPars.stimfilter, 'WNSTA')
        X_frame = (spfilter') * stim;
    end
end
frames = size(X_frame,2);
bpf    = GLMPars.bins_per_frame;
bins = bpf * frames;
X_bin  = repmat(X_frame, [ bpf,1]); 
X_bin = reshape(X_bin, [1, bins]);
clear frames bpf bins X_frame
%%
















end