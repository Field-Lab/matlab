% hack XV means cross validation
% here we use inputstats in our evaluation.. will work better for
% non-linearities  and NSEM 

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
function [X_frame,X_bin] = prep_stimcelldependentGPXV(GLMType, GLMPars, stimulus, inputstats, center_coord,STA, SU_filter, timefilter, troubleshoot)

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
ROI_length      = GLMPars.stimfilter.ROI_length;
stimsize.width  = size(stimulus,1);
stimsize.height = size(stimulus,2);
stimsize.frames = size(stimulus,3);
ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);

fitmoviestats.mean     =  inputstats.mu_avgIperpix;
fitmoviestats.span     =  inputstats.range;
fitmoviestats.normmean =  inputstats.mu_avgIperpix / inputstats.range;


%first subunits!
if exist('SU_filter', 'var') && SU_filter(1)~=0
    for frame=1:stimsize.frames
        tempstim = double(stimulus(:,:,frame))/double(fitmoviestats.span) - double(fitmoviestats.normmean);
        tempstim=filter2(SU_filter,tempstim);
        stim(:,:,frame)=tempstim(ROIcoord.xvals, ROIcoord.yvals);
    end
else
    stim   = double(stimulus(ROIcoord.xvals, ROIcoord.yvals, :));
    stim   = stim / double(fitmoviestats.span);
    if strcmp(GLMType.nullpoint, 'mean')
        stim = stim - double(fitmoviestats.normmean);
    elseif strcmp(GLMType.nullpoint, 'min')
        stim = stim - min(stim(:));
    else
        error('you need to fill in how to account for stimulus with a different nullpoint')
    end
end

% convolve the time filter if that's whats happening
if exist('timefilter', 'var') && timefilter(1)~=0
    timefilter = reshape(timefilter, [1 1 30]);
    stim_temp = convn(stim, timefilter, 'full');
    stim = stim_temp(:,:,1:stimsize.frames);    
end

% THEN the SU nonlinearity
if exist('SU_filter', 'var') && SU_filter(1)~=0
    if strcmp(GLMType.Subunit_NL, 'exp')
        stim=exp(stim);
    elseif strcmp(GLMType.Subunit_NL, 'squared')
        stim=stim.^2;
    end
end

% Actually should just use this for SU in the future... -NB 
if isfield(GLMType, 'input_pt_nonlinearity') && GLMType.input_pt_nonlinearity
   % display('implementing nonlinearity')
    newstim = stim;
    if strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_aboutmean')
        par       = GLMPars.others.point_nonlinearity.increment_to_decrement;
        pos_mult  = (2*par) / (par + 1) ;
        neg_mult  =      2  / (par + 1) ;
        
        pos_stim          = find(stim > 0 );
        neg_stim          = find(stim < 0 );
        
        newstim(pos_stim) = pos_mult * (newstim(pos_stim)); 
        newstim(neg_stim) = neg_mult * (newstim(neg_stim));
    elseif strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_shiftmean')
        par = GLMPars.others.point_nonlinearity.increment_to_decrement;
        par_shift = GLMPars.others.point_nonlinearity.shiftmean; 

        pos_mult  = (2*par) / (par + 1) ;
        neg_mult  =      2  / (par + 1) ;
        
        pos_stim          = find(stim > par_shift);
        neg_stim          = find(stim < par_shift);
        
        newstim(pos_stim) = pos_mult * (newstim(pos_stim)-par_shift) + par_shift; 
        newstim(neg_stim) = neg_mult * (newstim(neg_stim)-par_shift) + par_shift;        
    elseif strcmp(GLMType.input_pt_nonlinearity_type, 'log')
        display('implenting log')
        newstim = stim + fitmoviestats.normmean + 1; % now back on 0 1 scale
        newstim = log(newstim);
        newstim = newstim - log(fitmoviestats.normmean+1);
    elseif strcmp(GLMType.input_pt_nonlinearity_type, 'exp')
        display('implenting exp')
        newstim = stim + fitmoviestats.normmean; % now back on 0 1 scale
        newstim = exp(newstim);
        newstim = newstim - exp(fitmoviestats.normmean); 
    elseif strcmp(GLMType.input_pt_nonlinearity_type,'piecelinear_fourpiece_eightlevels')
        newstim = stim + fitmoviestats.normmean;
        
        quartile_1 = find(newstim<=.25);
        quartile_2 = setdiff(find(newstim<=.5),quartile_1);
        
        quartile_4 = find(newstim>.75);
        quartile_3 = setdiff(find(newstim >.5),quartile_4);
        
        
        coeff = GLMPars.others.point_nonlinearity.coefficients;
        slope1 = coeff.slope_quartile_1;
        slope2 = coeff.slope_quartile_2;
        slope3 = coeff.slope_quartile_3;
        slope4 = coeff.slope_quartile_4;
        
        offset1 = 0;
        offset2 = .25* slope1;
        offset3 = .25*(slope1+slope2);
        offset4 = .25*(slope1+slope2+slope3);
        
        newstim(quartile_1) = slope1 * (newstim(quartile_1) - .00) + offset1;
        newstim(quartile_2) = slope2 * (newstim(quartile_2) - .25) + offset2;
        newstim(quartile_3) = slope3 * (newstim(quartile_3) - .50) + offset3;
        newstim(quartile_4) = slope4 * (newstim(quartile_4) - .75) + offset4;
        
        newstim = newstim - offset2;
        
    % New option added AKHeitman 205-07-14    
    elseif strcmp(GLMType.input_pt_nonlinearity_type, 'log_powerraise')
        newstim = stim + fitmoviestats.normmean; % 0 1 scale        
        %display('log_powerraise')
        coeff = 10^(GLMPars.others.point_nonlinearity.log_powerraise);
        newstim = newstim.^coeff;
        newstim = newstim - mean(newstim(:));
    else
        display('error, need to properly specifiy input non-linearity')
    end
    stim = newstim; clear newstim
end

%}
stim = reshape(stim, [ROI_length^2 , stimsize.frames]);

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
        eval(sprintf('print -dpdf %s/%s_prepstimcellGP_spfilterfromSTA.pdf', troubleshoot.plotdir, troubleshoot.name));
    end 
    spfilter  = spfilter';  % dimension [1,ROI_length]
end

clear stimsize ROIcoord ROI_length fitmoviestats 

%%
%%%%%%%%%%%%%%%%%%%%%%%
% Creat the final stim dependent input (X) to the optimization algorithm
% X_bin  ([rank, bins])
%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
    if strcmp(GLMPars.stimfilter.fixedSP_type, 'WNSTA')
        X_frame = (spfilter) * stim;
    end
elseif strcmp(GLMType.stimfilter_mode, 'rk1') || strcmp(GLMType.stimfilter_mode, 'rk2') || ...
        strcmp(GLMType.stimfilter_mode, 'rk2-ConductanceBased')||strcmp(GLMType.stimfilter_mode, 'rk1-newfit')||strcmp(GLMType.stimfilter_mode, 'fullrank')
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
