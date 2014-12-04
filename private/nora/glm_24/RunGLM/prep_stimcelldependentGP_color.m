% AKHeitman 2014-04-27
% GLMPars = GLMParams  or GLMPars = GLMParams(GLMType.specialchange_name)
% Sensitive to naming changes in GLMParams.
% Lon function call is   "spatialfilterfromSTA"
% Only saves stuff (and calls directories) if we are in a troubleshooting mode
% Heavily GLMType dependent computations will be carried out here
% Outsourcable computations will be made into their own functions
%troubleshoot optional
% need troublshoot.doit (true or false), 
%troubleshoot.plotdir,
% troubleshoot.name
function [X_frame,X_bin] = prep_stimcelldependentGP_color(GLMType, GLMPars, stimulus, center_coord,STA,troubleshoot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up GLMParams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cut the Movie down to ROI
% Normalize Movie and set nullpoint
% output of this section is stim
% stim in [xy, time] coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stimsize.width  = size(stimulus,1);
stimsize.height = size(stimulus,2);
stimsize.frames = size(stimulus,4);
ROI_length      = GLMPars.stimfilter.ROI_length;
ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);

%{
%first subunits!
if GLMType.Subunits
	for frame=1:stimsize.frames
		tempstim=conv2(double(stimulus(:,:,frame)),double(GLMPars.subunits),'same');
		stimulus(:,:,frame)=tempstim;
	end
end
%}

stim            = stimulus(ROIcoord.yvals, ROIcoord.xvals, :,:);

fitmoviestats.mean     =  double(mean(stimulus(:)));
fitmoviestats.minval   =  double(min(stimulus(:)));
fitmoviestats.maxval   =  double(max(stimulus(:)));
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
%{
if isfield(GLMType, 'input_pt_nonlinearity') && GLMType.input_pt_nonlinearity
   % display('implementing nonlinearity')
    newstim = stim;
    if strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_aboutmean')
        par       = GLMPars.others.point_nonlinearity.increment_to_decrement;
        pos_mult  = (2*par) / (par + 1) ;
        neg_mult  =      2  / (par + 1) ; 
        pos_stim          = find(stim > 0 );
        newstim(pos_stim) = pos_mult * (newstim(pos_stim)); 
        neg_stim          = find(stim < 0 );
        newstim(neg_stim) = neg_mult * (newstim(neg_stim));
    end
    
    if strcmp(GLMType.input_pt_nonlinearity_type, 'raisepower_meanafter')
        newstim = newstim - min(newstim(:));
        newstim = newstim / max(newstim(:));  % back to a 0 1 scale
        
        newstim = newstim .^ (GLMPars.others.point_nonlinearity.scalar_raisedpower);
        newstim = newstim - min(newstim(:));
        newstim = newstim / max(newstim(:));
        %newstim = newstim - mean(newstim(:));
    end
    
    if strcmp(GLMType.input_pt_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
        par       = GLMPars.others.point_nonlinearity.scalar_raisedpower_aboutmean;

        pos_stim          = find(stim > 0 );
        newstim(pos_stim) = ((newstim(pos_stim)) .*par); 
        neg_stim          = find(stim < 0 );
        newstim(neg_stim) = -( (abs(newstim(pos_stim))) .*par);
    end
    stim = newstim; clear newstim
end
%}
%{
if exist('troubleshoot','var') && troubleshoot.doit
    clf
    totalframes = stimsize.frames;
    frames = min(1000, totalframes);
    stimulus2   = stimulus(:,:,1:frames);
    stim2       = stim(:,:,1:frames);
    subplot(3,2,[3 5]); hist(double(stimulus2(:)),20); set(gca, 'fontsize', 10); title('histogram of full raw stim')
    subplot(3,2,[4 6]); hist(           stim2(:) ,20); set(gca, 'fontsize', 10); title('histogram of normed stim over ROI') 
    subplot(3,2,[1 2]);   set(gca, 'fontsize', 10); axis off
    c = 0;
    c=c+1; text(-.1, 1-0.1*c,sprintf('Trouble shooting: %s',troubleshoot.name));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Specifically: Stimulus Normalizing Component of "glm-nospace"' ));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Plot Date %s',datestr(clock)));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')) );
    
    orient landscape
    eval(sprintf('print -dpdf %s/%s_prepstimcellGP_stimnorm.pdf', troubleshoot.plotdir, troubleshoot.name));
end 
%}
stim = reshape(stim, [ROI_length^2 , 3, stimsize.frames]);

clear stimsize ROIcoord ROI_length fitmoviestats 

%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create spfilter from the STA 
% output: spfilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
    ROI_length      = GLMPars.stimfilter.ROI_length;
    stimsize.width  = size(stimulus,1);
    stimsize.height = size(stimulus,2); 
    ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
    for i_filter=1:3
        spfilter{i_filter} = spatialfilterfromSTA(squeeze(STA(:,:,i_filter,:)),ROIcoord.yvals,ROIcoord.xvals)';
    end
    %{
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
        eval(sprintf('print -dpdf %s/%s_prepstimcellGP_spfilterfromSTA.pdf', troubleshoot.plotdir, troubleshoot.name));
    end 
    
    %}
    
     % dimension [1,ROI_length]
end

clear stimsize ROIcoord ROI_length fitmoviestats 

%%
%%%%%%%%%%%%%%%%%%%%%%%
% Creat the final stim dependent input (X) to the optimization algorithm
% X_bin  ([rank, bins])
%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
    if strcmp(GLMPars.stimfilter.fixedSP_type, 'WNSTA')
        X_frame=0;
        for i_filter=1:3
            X_frame = X_frame+(spfilter{i_filter}) * squeeze(stim(:,i_filter,:));
        end
    end
elseif strcmp(GLMType.stimfilter_mode, 'rk1') || strcmp(GLMType.stimfilter_mode, 'rk2') || strcmp(GLMType.stimfilter_mode, 'fullrank')
    X_frame = stim;
else
    error('you need to tell prep_stimcelldependentGP how to process stim for your spatial filter')
end

frames = size(X_frame,2);
dim    = size(X_frame,1);

bpf    = GLMPars.bins_per_frame;
bins   = bpf * frames;
X_bin  = repmat(X_frame, [ bpf,1]); 
X_bin  = reshape(X_bin, [dim, bins]);

clear frames bpf bins 


end
