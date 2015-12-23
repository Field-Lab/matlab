function com = ei_com_(ei,positions,nlPoints,varargin)
% soma_from_ei_     identify soma location based on EI location
%
% usage:  com = ei_com_(ei,positions,nlPoints,<params>)
%
% arguments:       ei - ExF matrix, standard EI format, E = # electrodes, F = # frames
%           positions - Ex2 matrix, list of electrode positions (in microns)
%            nlPoints - number of points before the spike, usually in datarun.ei.nlPoints
%            <params> - struct or list of optional arguments (see below)
%
% outputs:     com - 2-element cell array
%                       com{1} - Fx2 matrix, negative center of mass at each frame
%                       com{2} - Fx2 matrix, positive center of mass at each frame
%                       
%
% optional params, their default values, and what they specify:
%
% frames           	[]             	frames in which to find COM, relative to the spike time.  e.g. +2 = 2 frames post spike.
%                                       if empty, find COM in all frames
% foa               []             	figure or axes to plot in. if 0, make new figure.
%                                     	if empty, don't plot.  if -1, plot in current.
% roi               []              region of interest in which to compute the COM
%                                       {'peak',radius}     - radius (in microns) around the peak electrode
%                                       {'list',electrodes} - list of electrodes to use
%                                       if empty, use all electrodes
% axon              []              coordinates of an axon to plot--why not?
%
%
% 2010-03  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('foa', []);
p.addParamValue('frames', []);
p.addParamValue('roi', [], @iscell);
p.addParamValue('axon', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% set up plot axes
plot_axes = set_up_fig_or_axes(params.foa);


% threshold for significant electrodes, in EI sigmas
sig_elec_thresh = 5;


% identify ROI
if ~isempty(params.roi)
    switch params.roi{1}
        case 'peak'
            % identify electrode with largest amplitude signal
            [junk,peak] = max(max(ei,[],2));
            % get distance of each electrode from this one
            dists = sqrt(sum((positions - repmat(positions(peak,:),size(positions,1),1)).^2,2));
            % only use electrodes within the radius
            roi = dists < params.roi{2};
            
        case 'list'
            % if a list of electrodes, set just those electrodes to true
            roi = false(size(ei,1),1);
            roi(params.roi{2}) = true;
    end
    
else
    % if no ROI is specified, use all electrodes
    roi = true(size(ei,1),1);
end


% choose frames
if isempty(params.frames)
    params.frames = 1:size(ei,2);
else
    params.frames = params.frames + nlPoints;
end


% initialize
com = cell(2,1);
com{1} = zeros(length(params.frames),2);
com{2} = zeros(length(params.frames),2);



% for each frame, compute center of mass of large electrodes
for ff=1:length(params.frames)

    % note which frame this is
    frame = params.frames(ff);
    
    % compute threshold
    thresh = sig_elec_thresh*robust_std(ei(:,frame));
    
    % compute NEGATIVE COM
    
    % identify significant electrodes in the ROI
    sig_elec_neg = ei(:,frame) < -thresh & roi;
    % compute the amplitude of the EI on significant electrodes
    ei_amp = abs(ei(sig_elec_neg,frame));
    % compute the center of mass
    com{1}(ff,:) = [...
        ei_amp'*positions(sig_elec_neg,1)/sum(ei_amp)...
        ei_amp'*positions(sig_elec_neg,2)/sum(ei_amp)];
    thresh = sig_elec_thresh*robust_std(ei(:,frame));
    
    % compute POSITIVE COM
    
    % identify significant electrodes in the ROI
    sig_elec_pos = ei(:,frame) > thresh & roi;
    % compute the amplitude of the EI on significant electrodes
    ei_amp = abs(ei(sig_elec_pos,frame));
    % compute the center of mass
    com{2}(ff,:) = [...
        ei_amp'*positions(sig_elec_pos,1)/sum(ei_amp)...
        ei_amp'*positions(sig_elec_pos,2)/sum(ei_amp)];
    
    % plot
    if ~isempty(plot_axes)
        cla
        %plot_ei_(ei(:,frame),positions,1,'cutoff',-1,'scale',1,'pos_color','r','neg_color','b')
        plot_ei_(ei,positions,frame,'cutoff',-1,'scale',1,'pos_color','r','neg_color','b')
        
        plot_ei_(sig_elec_neg~=0,positions, 0,'alpha',0,'pretty_axes',0,'max_scale',1,'pos_color','c','neg_color','c','scale',1,'cutoff',-1)
        plot_ei_(sig_elec_pos~=0,positions, 0,'alpha',0,'pretty_axes',0,'max_scale',1,'pos_color','m','neg_color','m','scale',1,'cutoff',-1)
        plot(com{1}(1:ff,1),com{1}(1:ff,2),'.-','color',[.5 .5 1])
        plot(com{2}(1:ff,1),com{2}(1:ff,2),'.-','color',[1 .5 .5])
        title(num2str(frame))
        
        % add axon
        if ~isempty(params.axon)
            points = traced_cell_points(params.axon(1:2,:),params.axon(2:end,:));
            plot(points(:,1),points(:,2),'g')
        end
        
        %pause
        drawnow
    end
end



% post processing to find soma
if 0
    % compute how much the COM moved between each frame
    temp=diff(com);
    delta_com = sqrt(temp(:,1).^2 + temp(:,2).^2);


    % identify the longest contiguous set of frames where the com moved < com_jitter_thresh

    % identify jumps < com_jitter_thresh
    short_jumps = delta_com < com_jitter_thresh;

    % if there are none, return empty
    if isempty(short_jumps)
        soma_center = [];
        return
    end
    
    % ...

end

