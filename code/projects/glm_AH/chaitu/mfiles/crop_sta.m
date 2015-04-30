function [stas rstas cpidces box_m box_n ACmat] = crop_sta(sp_times,X,basepars,stimpars,rotateFlag,trialduration)

% Compute the STA for this cell using ALL spikes
N = basepars.Nneurons;

stas = cell(N,1);%zeros(size(X,1)*basepars.Mk,N);
rstas = cell(N,1);%zeros(basepars.n*basepars.Mk,N);
cpidces = cell(N,1);
1;

X = X - mean(X(:)); % subtract the mean of the stimulus before computing STA


if (isfield(stimpars,'trialtime'))
    trialduration = stimpars.trialtime;
end

if (~exist('trialduration','var'))
    trialduration = -1;
end
1;
for j=1:N
    
    if (iscell(sp_times))
        fprintf('Computing STA for cell %d/%d using %d spikes....',j,N,length(sp_times{j}));
        if (size(X,3)>1)
            
            sta = zeros(size(X,1),basepars.Mk);
            
            for k=1:size(X,3)
                sta = sta + 1/k.*fast_sta(X(:,:,j),get_spike_segment(sp_times{j},(k-1)*stimpars.trialtime,k*stimpars.trialtime,'rel'),basepars.Mk,stimpars.dt,0);
            end
        elseif (trialduration > 0)
            sta = fast_sta_rep(X,sp_times{j},trialduration,basepars.Mk,stimpars.dt,0);
        else
            sta = fast_sta(X,sp_times{j},basepars.Mk,stimpars.dt,0);
        end
    else
        fprintf('Computing STA for cell %d/%d using %d spikes....',j,N,length(sp_times(:,j)));
        
        if (size(X,3)>1)
            sta = zeros(size(X,1),basepars.Mk);
            for k=1:size(X,3)
                sta = sta + 1/k.*fast_sta(X(:,:,j),get_spike_segment(sp_times(:,j),(k-1)*stimpars.trialtime,k*stimpars.trialtime,'rel'),basepars.Mk,stimpars.dt,0);
            end
        elseif (trialduration > 0)
            sta = fast_sta_rep(X,sp_times(:,j),trialduration,basepars.Mk,stimpars.dt,0);
        else
            sta = fast_sta(X,sp_times(:,j),basepars.Mk,stimpars.dt,0);
        end
    end
    stas{j} = sta(:); % do NOT rotate out correlations before cropping
    fprintf('done.\n');
    
    1;
    
    if 0 % Play the STA if you want
        sta2play = (sta);
        sta2play = sta2play - min(sta2play(:));
        sta2play = sta2play./max(abs(sta2play(:)));
        play_sta(sta2play,basepars.stim_height,basepars.stim_width,basepars.Mk);
    end
    1;
    % Find cropping indices determined by the STAs and basepars
    %zmeanSTA = reshape(sum((sta - repmat(mean(sta,2),1,size(sta,2))).^2,2),basepars.stim_height,basepars.stim_width);
    zmeanSTA = reshape(sum((sta - repmat(mean(sta,2),1,size(sta,2))).^2,2),basepars.stim_width,basepars.stim_height);
    if (basepars.n > 1)
        if (isfield(basepars,'nocropflag') && basepars.nocropflag)
            cpidx = (1:basepars.n)'; box_m = sqrt(basepars.n); box_n = sqrt(basepars.n); % HACKKKK CHANGGGEGEEE
        else
            [cpidx box_m box_n] = get_crop_idx(zmeanSTA, sqrt(basepars.n),1); % get the relevant indices
        end
    else
        cpidx = [1];
        box_m = 1;
        box_n = 1;
    end
    1;
    if (isempty(cpidx)) % this is bad
        rstas{j} = 0;
        p0 = 0;
        ACmat = 0;
        cpidces = [];
        return;
    end

    cpidces{j} = cpidx;
    if 0
        sta2play = sta(cpidx,:);
        sta2play = sta2play - min(sta2play(:));
        sta2play = sta2play./max(abs(sta2play(:)));
        play_sta(sta2play,box_m,box_n,basepars.Mk);
    end
    1;
    sta_cropped = sta(cpidx,:);

    if rotateFlag
        %Rotate out correlations
        fprintf('Computing autocorrelation matrix\n');
        [rsta ACmat]= decorr_sta(X(cpidx,:),sta_cropped);
    else
        rsta = sta_cropped;
        rstas{j} = rsta(:);
        ACmat = eye(numel(rsta));
    end
    
end
