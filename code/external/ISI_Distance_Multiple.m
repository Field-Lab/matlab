% Calculates ISI-Distance from multiple spike trains "spikes" (matrix with different spike trains as rows)
% (new improved version that assigns just one value for each interval in the pooled spike train and is computationally faster and more memory-efficient)
%
% Information can be found online under "http://inls.ucsd.edu/~kreuz/Source-Code/Spike-Sync.html" and/or in
% Kreuz T, Haas JS, Morelli A, Abarbanel HDI, Politi A: Measuring spike train synchrony.
% J Neurosci Methods 165, 151 (2007)
% and
% Kreuz T, Chicharro D, Andrzejak RG, Haas JS, Abarbanel HDI: Measuring
% multiple spike train synchrony.
% J Neurosci Methods (submitted, 2009).
%
% Example call:
%
% num_trains=10; num_spikes=100;
% spikes=zeros(num_trains,num_spikes);
% spikes(1,1:num_spikes)=[1:num_spikes];          % first spike train (completely periodic)
% for trc=2:num_trains
%      spikes(trc,1:num_spikes-trc)=spikes(1,1:num_spikes-trc)+(rand(1,num_spikes-trc)-0.5);    % further spike trains (first spike train with some random jitter)
% end
% tmin=0; tmax=num_spikes+1; dt=0.001;            % (recording) edges and sampling interval
% id=ISI_Distance_Multiple(spikes,tmin,tmax,dt)   % id gives you the mean ISI, the mean rate, the averaged bivariate ISI-Distance and the ISI-Diversity 
%
% For questions and comment_strings please contact me at "tkreuz (at) ucsd.edu".

function isi_dist_multi=ISI_Distance_Multiple(spikes,tmin,tmax,dt)

max_mat_size=10000000;

if nargin<4
    % dt=0.001;           % Define your standard value (has to be smaller than the minimum ISI, ISIs should be multiples of dt)
    dt=f_get_dt(spikes);  % or get it automatically from the precision of the spikes
end

result_mode=15;           % +1:Mean ISI,+2:Mean rate,+4:Averaged bivariate ISI-Distance,+8:ISI-Diversity
weight_mode=1;            % 1-time weight_mode (paper-version), 2-spike weight_mode   (used for both calculation of the ISI-distance as well as the moving average)
edge_mode=1;              % 1-fixed interval,2:outer spikes,3-inner spikes
log_mode=0;               % 0-no log-transform,1-log-transform

window_mode=1;            % 1-all (recording limits), 2-smaller analysis window
if window_mode==2                              % if window_mode=2 select limits of the analysis window (wmin and wmax) here
    wmin=tmin+0.2*(tmax-tmin);
    wmax=tmax-0.2*(tmax-tmin);
    wmin=floor(wmin/dt)*dt;
    wmax=floor(wmax/dt)*dt;
    if wmin>=wmax
        disp(' '); disp(' '); error ('********** Error: W_max must be larger than W_min !!!!!')
    end
    if wmin<tmin
        disp(' '); disp(' '); error ('********** Error: W_min must not be smaller than T_min !!!!!')
    end
    if wmax>tmax
        disp(' '); disp(' '); error ('********** Error: W_max must not be larger than T_max !!!!!')
    end
end

plot_mode=1;             % plotting figure? (0-no,1-yes)
colors='mkbrg';          % Colors for 1-Stimulus, 2-Spike Trains / ISI, 3-Mean values, 4-Moving averages of mean values, 5-Window
num_fig=1;               % number of (first) figure
title_string='Multiple-Example :';     % will be part of the title
comment_string='';       % will be part of the title and the filename
timeunit_string='[ms]';  % time unit
patch_mode=0;            % patching? (0-no,1-yes)

print_mode=1;            % saving to postscript file? (0-no,1-yes)
filename='Multiple-ISI-Distance-Test';  % filename of this postscript file


intmean_mao=10;     % order of the moving average for the intervals
ratemean_mao=10;   % order of the moving average for the rates
pa_mao=10;         % order of the moving average for the pairwise average
cv_mao=10;         % order of the moving average for the coefficent of variation

% Choose order and size of the subplots; for multiple figures use a matrix ( in which each row corresponds to one figure)
% Use "subplot_posi_mat" to select the order of the subplots. Use 0 if a subplot is needed.
% Make sure that the range from 1 to the desired number of subplots is covered
% Use "subplot_size_mat" to select the size of the selected subplots.
% Its order refers to the subplots selected in subplot_posi_mat, not to the reference-order
% (i.e., its length should be equal to the number of subplot selected in 'subplot_posi_mat')

%
% The subplot "Stimulus" allows you to add any representation of the stimulus that triggered the spike trains (put 0 if not needed)

% Reference-Order (for position variable only): Stimulus Spikes ISI ISI-Mean ISI-Mean-MA Rate Rate-Mean Rate-Mean-MA A A-MA CV CV-MA
% only 3, 4, 5 and 6, 7, 8 and 9, 10 as well as 11, 12 use the same scale and thus can be the same number

%subplot_posi_mat=[0 1 2 2 2 3 3 3 4 4 5 5];
%subplot_size_mat=[1 1 1 1 1];

subplot_posi_mat=[0 1 0 2 0 0 3 0 4 0 5 0];
subplot_size_mat=[1 1 1 1 1];

% This example (if uncommented) shows all plots (except for the stimulus) separated
% subplot_posi_mat=[0 1 2 3 4 5 6 7 8 9 10 11];    % Order of suplots
% subplot_size_mat=[1 1 1 1 1 1 1 1 1 1 1];        % Relative size of these subplots
% This example (if uncommented) shows all plots (except for the stimulus) with maximum overlap
% subplot_posi_mat=[0 1 2 2 2 3 3 3 4 4 5 5];
% subplot_size_mat=[1 1 1 1 1];

num_trains=size(spikes,1);
num_pairs=(num_trains*num_trains-num_trains)/2;
num_pspikes=zeros(1,num_trains);
for trac=1:num_trains
    if any(spikes(trac,:))
        num_pspikes(trac)=find(spikes(trac,:),1,'last');
    end
end
max_num_pspikes=max(num_pspikes);

pspikes=zeros(num_trains,max_num_pspikes);          % original spikes used for plotting
for trac=1:num_trains
    pspikes(trac,1:num_pspikes(trac))=floor(sort(spikes(trac,1:num_pspikes(trac)))/dt)*dt;
end

if nargin<3 && exist('tmin','var')==0  && exist('tmax','var')==0
    tmin=min(min(pspikes));
    tmax=max(max(pspikes));
    trange=tmax-tmin;
    tmin=floor((tmin-0.02*trange)/dt)*dt;
    tmax=floor((tmax+0.02*trange)/dt)*dt;
else
    tmin=floor(tmin/dt)*dt;
    tmax=floor(tmax/dt)*dt;
end

if edge_mode==1                                                               % fixed interval
    itmin=tmin;
    itmax=tmax;
elseif edge_mode==2                                                           % outer spikes (overall)
    itmin=min(min(min(pspikes(pspikes~=0))));
    itmax=max(max(max(pspikes(pspikes~=0))));
elseif edge_mode==3                                                           % inner spikes (overall)
    itmin=max(pspikes(1:num_trains,1));
    itmax=inf;
    for trac=1:num_trains
        if num_pspikes(trac)>0
            if pspikes(trac,num_pspikes(trac))<itmax
                itmax=pspikes(trac,num_pspikes(trac));
            end
        end
    end
end
itmin=floor(itmin/dt)*dt;
itmax=floor(itmax/dt)*dt;
itrange=itmax-itmin;

if window_mode==2
    wmin=max([itmin wmin]);
    wmax=min([itmax wmax]);
end

num_ispikes=zeros(1,num_trains);
ispikes=zeros(num_trains,max_num_pspikes+2);
for trac=1:num_trains
    dspikes=pspikes(trac,1:num_pspikes(trac));
    dspikes=dspikes(dspikes>=itmin & dspikes<=itmax);
    if any(dspikes)
        num_ispikes(trac)=find(dspikes,1,'last')+2;
        if isempty(find(dspikes(1:num_ispikes(trac)-2)==itmin,1))
            ispikes(trac,1:num_ispikes(trac))=[itmin dspikes(1:num_ispikes(trac)-2) itmax];
        else
            dummy=dspikes(1:num_ispikes(trac)-2);
            dummy2=[itmin dummy(dummy~=itmin) itmax];
            num_ispikes(trac)=length(dummy2);
            ispikes(trac,1:num_ispikes(trac))=dummy2;
        end
    else
        num_ispikes(trac)=2;
        ispikes(trac,1:num_ispikes(trac))=[itmin itmax];
    end
end
max_num_ispikes=max(num_ispikes);
ispikes=ispikes(1:num_trains,1:max_num_ispikes);

num_coins=zeros(num_trains,max_num_ispikes);
num_uspikes=zeros(1,num_trains);
uspikes=zeros(num_trains,max_num_ispikes);
for trac=1:num_trains
    uispikes=unique(ispikes(trac,1:num_ispikes(trac)));
    for uic=1:length(uispikes)-1
        num_coins(trac,uic)=sum(ispikes(trac,1:num_ispikes(trac))==uispikes(uic));
    end
    num_uspikes(trac)=length(uispikes);
    uspikes(trac,1:num_uspikes(trac))=uispikes;
end
max_num_uspikes=max(num_uspikes);
uspikes=uspikes(1:num_trains,1:max_num_uspikes);

for trac=1:num_trains
    uspikes(trac,num_uspikes(trac)+1:max_num_uspikes)=987654321.123456789;
end

[all_spikes,all_indy]=sort(reshape(uspikes',1,numel(uspikes)));
all_labels=fix((all_indy-1)/max_num_uspikes)+1;
all_labels=all_labels(all_spikes~=987654321.123456789);
all_spikes=all_spikes(all_spikes~=987654321.123456789);
all_labels(1:num_trains)=0;
all_labels(end-num_trains+1:end)=0;

all_spikes2=all_spikes(num_trains:end-num_trains+1);
all_labels2=all_labels(num_trains:end-num_trains+1);
all_isi2=diff(all_spikes2);
num_all_isi2=length(all_isi2);
indy=find(all_spikes2==itmin,1,'last'):find(all_spikes2==itmax,1,'first')-1;

isis=zeros(num_trains,max_num_uspikes-1);
for trac=1:num_trains
    isis(trac,1:num_uspikes(trac)-1)=diff(uspikes(trac,1:num_uspikes(trac)));
    isis(trac,isis(trac,:)~=0)=isis(trac,isis(trac,:)~=0)./num_coins(trac,isis(trac,:)~=0);
end

ints=zeros(num_trains,num_all_isi2);
rates=zeros(num_trains,num_all_isi2);
for trac=1:num_trains
    ivs=[1 find(all_labels2==trac)];
    ive=[ivs(2:num_uspikes(trac)-1)-1 num_all_isi2];
    for ic=1:num_uspikes(trac)-1
        ints(trac,ivs(ic):ive(ic))=isis(trac,ic);
    end
    rates(trac,ints(trac,:)>0)=1./ints(trac,ints(trac,:)>0);
end

mean_ints=mean(ints);
mean_rates=mean(rates);
if weight_mode==1                                                          % time-weighted (paper-version)
    multi_mean_isi=sum(mean_ints.*all_isi2)/sum(all_isi2);
    multi_mean_rate=sum(mean_rates.*all_isi2)/sum(all_isi2);
else
    multi_mean_isi=mean(mean_ints);
    multi_mean_rate=mean(mean_rates);
end

isi_dist_multi=[];
if mod(result_mode,2)>0
    isi_dist_multi=[isi_dist_multi multi_mean_isi];
end
if mod(result_mode,4)>1
    isi_dist_multi=[isi_dist_multi multi_mean_rate];
end

if any(subplot_posi_mat(:,9:10))>0 || mod(result_mode,8)>3
    mat_size=num_pairs*num_all_isi2;
    if mat_size>max_mat_size
        if max_mat_size<num_pairs
            max_mat_size=num_pairs;
        end
        max_num_all_isi2=fix(max_mat_size/num_pairs);
        num_runs=ceil(num_all_isi2/max_num_all_isi2);
        nai_runs=[max_num_all_isi2*ones(1,num_runs-1) num_all_isi2-max_num_all_isi2*(num_runs-1)];
        run_ends=cumsum(nai_runs);
        run_starts=[1 run_ends(1:end-1)+1];
    else
        num_runs=1;
        nai_runs=num_all_isi2;
        run_starts=1;
        run_ends=num_all_isi2;
    end

    pa_ints=zeros(1,num_all_isi2);
    run_sum_all_isi2=zeros(1,num_runs);
    run_multi_pa_isi_dist=zeros(1,num_runs);
    for ruc=1:num_runs
        run_xy_isi_dist=zeros(num_pairs,nai_runs(ruc));
        run_range=run_starts(ruc):run_ends(ruc);
        pac=0;
        for trac=1:num_trains-1
            for trac2=trac+1:num_trains
                pac=pac+1;
                dum=find(ints(trac,run_starts(ruc):run_ends(ruc))<ints(trac2,run_starts(ruc):run_ends(ruc)));
                run_xy_isi_dist(pac,dum)=ints(trac,run_range(dum))./ints(trac2,run_range(dum))-1;
                dum2=find(ints(trac,run_starts(ruc):run_ends(ruc))>=ints(trac2,run_starts(ruc):run_ends(ruc)) & ints(trac,run_starts(ruc):run_ends(ruc))~=0);
                run_xy_isi_dist(pac,dum2)=-(ints(trac2,run_range(dum2))./ints(trac,run_range(dum2))-1);
            end
        end
        run_pa_ints=mean(abs(run_xy_isi_dist),1);
        run_indy=indy(indy>=run_starts(ruc) & indy<=run_ends(ruc))-run_starts(ruc)+1;
        run_sum_all_isi2(ruc)=sum(all_isi2(run_starts(ruc)-1+run_indy));
        run_multi_pa_isi_dist(ruc)=sum(abs(run_pa_ints(run_indy).*all_isi2(run_starts(ruc)-1+run_indy)))/run_sum_all_isi2(ruc); % mean over pairwise distances
        pa_ints(run_range)=run_pa_ints;
    end

    if weight_mode==1                                                          % time-weighted (paper-version)
        all_multi_pa_isi_dist=sum(pa_ints.*all_isi2)/sum(all_isi2);
        if window_mode==1
            multi_pa_isi_dist=all_multi_pa_isi_dist;
        else
            first_winspike_ind=find(all_spikes2>=wmin,1,'first');
            last_winspike_ind=find(all_spikes2<=wmax,1,'last');
            if first_winspike_ind>last_winspike_ind
                win_all_isi2=wmax-wmin;
                win_pa_ints=pa_ints(last_winspike_ind:first_winspike_ind-1);
            else
                win_all_isi2=all_isi2(first_winspike_ind:last_winspike_ind-1);
                win_pa_ints=pa_ints(first_winspike_ind:last_winspike_ind-1);
                if wmin<all_spikes2(first_winspike_ind);
                    win_all_isi2=[all_spikes2(first_winspike_ind)-wmin win_all_isi2];
                    win_pa_ints=[pa_ints(first_winspike_ind-1) win_pa_ints];
                end
                if wmax>all_spikes2(last_winspike_ind);
                    win_all_isi2=[win_all_isi2 wmax-all_spikes2(last_winspike_ind)];
                    win_pa_ints=[win_pa_ints pa_ints(last_winspike_ind)];
                end
            end
            win_multi_pa_isi_dist=sum(win_pa_ints.*win_all_isi2)/sum(win_all_isi2);
            multi_pa_isi_dist=win_multi_pa_isi_dist;
        end
    else                                                                       % spike-weighted
        all_multi_pa_isi_dist=mean(pa_ints);
        if window_mode==1
            multi_pa_isi_dist=all_multi_pa_isi_dist;
        else
            first_winspike_ind=find(all_spikes2>=wmin,1,'first');
            last_winspike_ind=find(all_spikes2<=wmax,1,'last');
            win_pa_ints=pa_ints(first_winspike_ind:last_winspike_ind-1);
            if wmin<all_spikes2(first_winspike_ind);
                win_pa_ints=[pa_ints(first_winspike_ind-1) win_pa_ints];
            end
            if wmax>all_spikes2(last_winspike_ind);
                win_pa_ints=[win_pa_ints pa_ints(last_winspike_ind)];
            end
            win_multi_pa_isi_dist=mean(win_pa_ints);
            multi_pa_isi_dist=win_multi_pa_isi_dist;
        end
    end

    result_string=['  D_I^a = ',num2str(multi_pa_isi_dist,3)];
    if mod(result_mode,8)>3
        isi_dist_multi=[isi_dist_multi multi_pa_isi_dist];
    end
else
    result_string='';
end

if any(subplot_posi_mat(:,11:12))>0 || mod(result_mode,8)>3
    if log_mode==0
        cv_ints=std(ints)./mean(ints);
    else
        log_ints=log(round(ints/dt));
        cv_ints=mean(log_ints).*(exp(std(log_ints)./mean(log_ints))-1);
    end
    if weight_mode==1                                                          % time-weighted (paper-version)
        all_multi_cv_isi_dist=sum(cv_ints(indy).*all_isi2)/sum(all_isi2);
        if window_mode==1
            multi_cv_isi_dist=all_multi_cv_isi_dist;
        else
            first_winspike_ind=find(all_spikes2>wmin,1,'first');
            last_winspike_ind=find(all_spikes2<wmax,1,'last');
            if first_winspike_ind>last_winspike_ind
                win_all_isi2=wmax-wmin;
                win_cv_ints=cv_ints(last_winspike_ind:first_winspike_ind-1);
            else
                win_all_isi2=all_isi2(first_winspike_ind:last_winspike_ind-1);
                win_cv_ints=cv_ints(first_winspike_ind:last_winspike_ind-1);
                if wmin<all_spikes2(first_winspike_ind);
                    win_all_isi2=[all_spikes2(first_winspike_ind)-wmin win_all_isi2];
                    win_cv_ints=[cv_ints(first_winspike_ind-1) win_cv_ints];
                end
                if wmax>all_spikes2(last_winspike_ind);
                    win_all_isi2=[win_all_isi2 wmax-all_spikes2(last_winspike_ind)];
                    win_cv_ints=[win_cv_ints cv_ints(last_winspike_ind)];
                end
            end
            win_multi_cv_isi_dist=sum(win_cv_ints.*win_all_isi2)/sum(win_all_isi2);
            multi_cv_isi_dist=win_multi_cv_isi_dist;
        end
    else                                                                       % spike-weighted
        all_multi_cv_isi_dist=mean(cv_ints(indy));
        if window_mode==1
            multi_cv_isi_dist=all_multi_cv_isi_dist;
        else
            first_winspike_ind=find(all_spikes2>wmin,1,'first');
            last_winspike_ind=find(all_spikes2<wmax,1,'last');
            win_cv_ints=cv_ints(first_winspike_ind:last_winspike_ind-1);
            if wmin<all_spikes2(first_winspike_ind);
                win_cv_ints=[cv_ints(first_winspike_ind-1) win_cv_ints];
            end
            if wmax>all_spikes2(last_winspike_ind);
                win_cv_ints=[win_cv_ints cv_ints(last_winspike_ind)];
            end
            win_multi_cv_isi_dist=mean(win_cv_ints);
            multi_cv_isi_dist=win_multi_cv_isi_dist;
        end
    end

    if any(subplot_posi_mat(:,9:10))>0
        result_string=[result_string,'  ;  D_I^m = ',num2str(multi_cv_isi_dist,3)];
    else
        result_string=['  D_I^m = ',num2str(multi_cv_isi_dist,3)];
    end
    if mod(result_mode,16)>7
        isi_dist_multi=[isi_dist_multi multi_cv_isi_dist];
    end
end

if plot_mode>0
    if edge_mode==1                                                              % limits for plotting
        pmin=itmin-0.02*itrange;
        pmax=itmax+0.02*itrange;
    else
        pmax2=-inf;
        for trac=1:num_trains
            if num_pspikes(trac)>0
                pmax2=max([pmax2 pspikes(trac,num_pspikes(trac))]);
            end
        end
        pmin2=min(pspikes(1:num_trains,1));
        prange2=pmax2-pmin2;
        pmin=pmin2-0.02*prange2;
        pmax=pmax2+0.02*prange2;
    end
    maxintval=max(max(ints(:,indy)));
    maxrateval=max(max(rates(:,indy)));


    pints_support=zeros(num_trains,max_num_uspikes-2);
    pints=zeros(num_trains,max_num_uspikes-2);
    prates=zeros(num_trains,max_num_uspikes-2);
    for trac=1:num_trains
        pints_support(trac,1:2*num_uspikes(trac)-2)=sort([uspikes(trac,1:num_uspikes(trac)) uspikes(trac,2:num_uspikes(trac)-1)-dt]);
        pints(trac,1:2*num_uspikes(trac)-2)=reshape([diff(uspikes(trac,1:num_uspikes(trac))); diff(uspikes(trac,1:num_uspikes(trac)))],1,num_uspikes(trac)*2-2);
        prates(trac,1:2*num_uspikes(trac)-2)=1./pints(trac,1:2*num_uspikes(trac)-2);
    end
    pints(pints>maxintval)=maxintval;
    prates(prates>maxrateval)=maxrateval;
    pratios_support=sort([all_spikes2 all_spikes2(2:num_all_isi2)-dt]);

    if any(subplot_posi_mat(:,4))
        pmean_ints=reshape([mean_ints; mean_ints],1,num_all_isi2*2);
    end
    if any(subplot_posi_mat(:,5))
        if weight_mode==1
            mean_ints_ma=f_weighted_moving_average(mean_ints,all_isi2,intmean_mao);      % time-weighted (paper-version)
        else
            mean_ints_ma=f_moving_average(mean_ints,intmean_mao);                        % spike-weighted
        end
        pmean_ints_ma=reshape([mean_ints_ma; mean_ints_ma],1,num_all_isi2*2);
    end
    if any(subplot_posi_mat(:,7))
        pmean_rates=reshape([mean_rates; mean_rates],1,num_all_isi2*2);
    end
    if any(subplot_posi_mat(:,8))
        if weight_mode==1
            mean_rates_ma=f_weighted_moving_average(mean_rates,all_isi2,ratemean_mao);      % time-weighted (paper-version)
        else
            mean_rates_ma=f_moving_average(mean_rates,ratemean_mao);                        % spike-weighted
        end
        pmean_rates_ma=reshape([mean_rates_ma; mean_rates_ma],1,num_all_isi2*2);
    end
    if any(subplot_posi_mat(:,9)) || any(subplot_posi_mat(:,10))
        ppa_ints=reshape([pa_ints; pa_ints],1,num_all_isi2*2);
        if any(subplot_posi_mat(:,10))
            if weight_mode==1
                pa_ints_ma=f_weighted_moving_average(pa_ints,all_isi2,pa_mao);      % time-weighted (paper-version)
            else
                pa_ints_ma=f_moving_average(pa_ints,pa_mao);                        % spike-weighted
            end
            ppa_ints_ma=reshape([pa_ints_ma; pa_ints_ma],1,num_all_isi2*2);
        end
    end
    if any(subplot_posi_mat(:,11)) || any(subplot_posi_mat(:,12))
        pcv_ints=reshape([cv_ints; cv_ints],1,num_all_isi2*2);
        if any(subplot_posi_mat(:,12))
            if weight_mode==1
                cv_ints_ma=f_weighted_moving_average(cv_ints,all_isi2,cv_mao);      % time-weighted (paper-version)
            else
                cv_ints_ma=f_moving_average(cv_ints,cv_mao);                        % spike-weighted
            end
            pcv_ints_ma=reshape([cv_ints_ma; cv_ints_ma],1,num_all_isi2*2);
        end
    end

    num_figs=size(subplot_posi_mat,1);
    for plc=1:num_figs
        subplot_posi=subplot_posi_mat(plc,:);
        subplot_size=subplot_size_mat(plc,:);

        if ~any(subplot_posi==1)
            subplot_posi=(subplot_posi-min(subplot_posi(subplot_posi>0))+1).*(subplot_posi>0);
        end

        relsubplot_size=subplot_size(1:find(subplot_size>0,1,'last'));
        if length(relsubplot_size)~=max(subplot_posi)
            relsubplot_size=ones(1,max(subplot_posi));
        end

        singles=zeros(1,max(subplot_posi));
        for suc=1:max(subplot_posi)
            singles(suc)=find(subplot_posi==suc,1,'first');                               % new subplot
        end
        num_subplots=length(singles);

        doubles=intersect(setxor(1:length(subplot_posi),singles),find(subplot_posi>0));     % repeated subplot
        doublesref=singles(subplot_posi(doubles));                                                % corresponding new subplot

        regplotsize=num_subplots*1.1;
        normsubplot_size=relsubplot_size/sum(relsubplot_size);
        dumsubplot_size=normsubplot_size*regplotsize;
        subplot_start2=cumsum(dumsubplot_size);
        subplot_size2=diff([0 subplot_start2]);

        subplot_size=zeros(1,length(subplot_posi));               % normalized size of subplots
        subplot_size(singles)=subplot_size2;
        subplot_size(doubles)=subplot_size(doublesref);

        subplot_start=zeros(1,length(subplot_posi));               % normalized position of subplots
        subplot_start(singles)=subplot_start2;
        subplot_start(doubles)=subplot_start(doublesref);

        figure(num_fig-1+plc); clf; hold on;
        set(gcf,'Position',[5 39 1912 1072])
        set(gcf,'Name',['ISI-Multi-Distance : ',title_string,'   ',comment_string])
        xlim([pmin pmax])
        ylim([0 regplotsize])
        xl=xlim; yl=ylim;
        for spc=1:length(subplot_start)
            if subplot_posi(spc)>0 && subplot_start(spc)~=regplotsize && subplot_start(spc)~=0
                line(xl,(yl(2)-subplot_start(spc))*ones(1,2),'Color','k','LineStyle','-','LineWidth',2)
            end
        end
        yt=[]; ytl=[];

        if subplot_posi(1)>0                                                                          % Stimulus
            % here you can plot the stimulus. The example below (if uncomment_stringed) shows a sine wave.
            % stim=0.5+sin((0:1/((itmax-itmin)/dt):1)*2*pi)/2;
            % plot(itmin:dt:itmax,yl(2)-subplot_start(1)+0.05+stim/1.1*subplot_size(1),['.-',colors(1)])
            % max_stim_val=1;
            % stim_lab=[-1 0 1];
            % yt=[yt yl(2)-subplot_start(1)+(0.05+[0 0.5 1]/max_stim_val)/1.1*subplot_size(1)];
            % ytl=[ytl stim_lab];
            % text(xl(1)-0.095*(xl(2)-xl(1)),yl(2)-subplot_start(1)+0.75/1.1*subplot_size(1),'Stimulus','Color','k','FontSize',11,'FontWeight','bold')
            % line(xl,yl(2)-subplot_start(1)+0.05/1.1*subplot_size(1)*ones(1,2),'Color','k','LineStyle',':')
            % line(xl,yl(2)-subplot_start(1)+1.05/1.1*subplot_size(1)*ones(1,2),'Color','k','LineStyle',':')
        end

        if subplot_posi(2)>0                                                                          % Spikes
            for trac=1:num_trains
                for sc=1:num_pspikes(trac)
                    line(pspikes(trac,sc)*ones(1,2),yl(2)-subplot_start(2)+(0.05+(num_trains-1-(trac-1)+[0.05 0.95])/num_trains)/1.1*subplot_size(2),'Color',colors(2))
                end
            end
            spikesvals=[0 0.5 1]; num_vals=2;
            spikeslab=fliplr(f_lab(spikesvals*num_trains,num_vals,1,1));
            yt=[yt yl(2)-subplot_start(2)+(1.05-[spikeslab 0]/num_trains)/1.1*subplot_size(2)];
            ytl=[ytl spikeslab 0];
            text(xl(1)-0.09*(xl(2)-xl(1)),yl(2)-subplot_start(2)+0.75/1.1*subplot_size(2),'Spike','Color','k','FontSize',11,'FontWeight','bold')
            text(xl(1)-0.09*(xl(2)-xl(1)),yl(2)-subplot_start(2)+0.35/1.1*subplot_size(2),'trains','Color','k','FontSize',11,'FontWeight','bold')
            line(xl,yl(2)-subplot_start(2)+0.05/1.1*subplot_size(2)*ones(1,2),'Color','k','LineStyle',':')
            line(xl,yl(2)-subplot_start(2)+1.05/1.1*subplot_size(2)*ones(1,2),'Color','k','LineStyle',':')
        end

        if subplot_posi(3)>0                                                                          % ISI
            for trac=1:num_trains
                plot(pints_support(trac,1:2*num_uspikes(trac)-2),yl(2)-subplot_start(3)+(0.05+pints(trac,1:2*num_uspikes(trac)-2)/maxintval)/1.1*subplot_size(3),['-',colors(2)])
            end
            intvals=[0.5 1]; num_vals=2;
            intlab=f_lab(intvals*maxintval,num_vals,1,1);
            yt=[yt yl(2)-subplot_start(3)+(0.05+[0 intlab]/maxintval)/1.1*subplot_size(3)];
            ytl=[ytl 0 intlab];
            text(xl(1)-0.08*(xl(2)-xl(1)),yl(2)-subplot_start(3)+0.75/1.1*subplot_size(3),'ISI','Color','k','FontSize',11,'FontWeight','bold')
            line(xl,yl(2)-subplot_start(3)+0.05/1.1*subplot_size(3)*ones(1,2),'Color','k','LineStyle',':')
            line(xl,yl(2)-subplot_start(3)+1.05/1.1*subplot_size(3)*ones(1,2),'Color','k','LineStyle',':')
        end

        if subplot_posi(4)>0                                                                          % ISI-Mean
            if subplot_posi(4)~=subplot_posi(3)
                maxmeanintval=max(pmean_ints);
                intvals=[0.5 1]; num_vals=2;
                intlab=f_lab(intvals*maxmeanintval,num_vals,1,1);
                yt=[yt yl(2)-subplot_start(4)+(0.05+[0 intlab]/maxmeanintval)/1.1*subplot_size(4)];
                ytl=[ytl 0 intlab];
                text(xl(1)-0.09*(xl(2)-xl(1)),yl(2)-subplot_start(4)+0.75/1.1*subplot_size(4),'<ISI>','Color','k','FontSize',11,'FontWeight','bold')
                line(xl,yl(2)-subplot_start(4)+0.05/1.1*subplot_size(4)*ones(1,2),'Color','k','LineStyle',':')
                line(xl,yl(2)-subplot_start(4)+1.05/1.1*subplot_size(4)*ones(1,2),'Color','k','LineStyle',':')
                lw=1;
            else
                maxmeanintval=maxintval;
                lw=2;
            end
            plot(pratios_support,yl(2)-subplot_start(4)+(0.05+pmean_ints/maxmeanintval)/1.1*subplot_size(4),['-',colors(3)],'LineWidth',lw)
            %line([itmin itmax],yl(2)-subplot_start(4)+(0.05+multi_mean_isi/maxmeanintval)/1.1*subplot_size(4)*ones(1,2),'Color',colors(3),'LineStyle','--')
        end

        if subplot_posi(5)>0                                                                          % ISI-Mean-MA
            if ~any(subplot_posi(3:4)==subplot_posi(5))
                maxmeanintmaval=max(pmean_ints_ma);
                intvals=[0.5 1]; num_vals=2;
                intlab=f_lab(intvals*maxmeanintmaval,num_vals,1,1);
                yt=[yt yl(2)-subplot_start(5)+(0.05+[0 intlab]/maxmeanintmaval)/1.1*subplot_size(5)];
                ytl=[ytl 0 intlab];
                text(xl(1)-0.09*(xl(2)-xl(1)),yl(2)-subplot_start(5)+0.75/1.1*subplot_size(5),'<ISI>^*','Color','k','FontSize',11,'FontWeight','bold')
                line(xl,yl(2)-subplot_start(5)+0.05/1.1*subplot_size(5)*ones(1,2),'Color','k','LineStyle',':')
                line(xl,yl(2)-subplot_start(5)+1.05/1.1*subplot_size(5)*ones(1,2),'Color','k','LineStyle',':')
                lw=1;
            else
                if subplot_posi(5)==subplot_posi(3)
                    maxmeanintmaval=maxintval;
                else
                    maxmeanintmaval=maxmeanintval;
                end
                lw=2.5;
            end
            plot(pratios_support,yl(2)-subplot_start(5)+(0.05+pmean_ints_ma/maxmeanintmaval)/1.1*subplot_size(5),['-',colors(4)],'LineWidth',lw)
        end

        if subplot_posi(6)>0                                                                          % Rate
            for trac=1:num_trains
                plot(pints_support(trac,1:2*num_uspikes(trac)-2),yl(2)-subplot_start(6)+(0.05+prates(trac,1:2*num_uspikes(trac)-2)/maxrateval)/1.1*subplot_size(6),['-',colors(2)])
            end
            intvals=[0.5 1]; num_vals=2;
            intlab=f_lab(intvals*maxrateval,num_vals,1,1);
            yt=[yt yl(2)-subplot_start(6)+(0.05+[0 intlab]/maxrateval)/1.1*subplot_size(6)];
            ytl=[ytl 0 intlab];
            text(xl(1)-0.09*(xl(2)-xl(1)),yl(2)-subplot_start(6)+0.75/1.1*subplot_size(6),'1/ISI','Color','k','FontSize',11,'FontWeight','bold')
            line(xl,yl(2)-subplot_start(6)+0.05/1.1*subplot_size(6)*ones(1,2),'Color','k','LineStyle',':')
            line(xl,yl(2)-subplot_start(6)+1.05/1.1*subplot_size(6)*ones(1,2),'Color','k','LineStyle',':')
        end

        if subplot_posi(7)>0                                                                          % Rate-Mean
            if subplot_posi(7)~=subplot_posi(6)
                maxmeanrateval=max(pmean_rates);
                ratevals=[0.5 1]; num_vals=2;
                ratelab=f_lab(ratevals*maxmeanrateval,num_vals,1,1);
                yt=[yt yl(2)-subplot_start(7)+(0.05+[0 ratelab]/maxmeanrateval)/1.1*subplot_size(7)];
                ytl=[ytl 0 ratelab];
                text(xl(1)-0.08*(xl(2)-xl(1)),yl(2)-subplot_start(7)+0.75/1.1*subplot_size(7),'R','Color','k','FontSize',11,'FontWeight','bold')
                line(xl,yl(2)-subplot_start(7)+0.05/1.1*subplot_size(7)*ones(1,2),'Color','k','LineStyle',':')
                line(xl,yl(2)-subplot_start(7)+1.05/1.1*subplot_size(7)*ones(1,2),'Color','k','LineStyle',':')
                lw=1;
            else
                maxmeanrateval=maxrateval;
                lw=2;
            end
            plot(pratios_support,yl(2)-subplot_start(7)+(0.05+pmean_rates/maxmeanrateval)/1.1*subplot_size(7),['-',colors(3)],'LineWidth',lw)
            %line([itmin itmax],yl(2)-subplot_start(7)+(0.05+multi_mean_rate/maxrateval)/1.1*subplot_size(7)*ones(1,2),'Color',colors(3),'LineStyle','--')
        end

        if subplot_posi(8)>0                                                                          % Rate-Mean-MA
            if ~any(subplot_posi(6:7)==subplot_posi(8))
                maxmeanratemaval=max(pmean_rates_ma);
                ratevals=[0.5 1]; num_vals=2;
                ratelab=f_lab(ratevals*maxmeanratemaval,num_vals,1,1);
                yt=[yt yl(2)-subplot_start(8)+(0.05+[0 ratelab]/maxmeanratemaval)/1.1*subplot_size(8)];
                ytl=[ytl 0 ratelab];
                text(xl(1)-0.08*(xl(2)-xl(1)),yl(2)-subplot_start(8)+0.75/1.1*subplot_size(8),'R^*','Color','k','FontSize',11,'FontWeight','bold')
                line(xl,yl(2)-subplot_start(8)+0.05/1.1*subplot_size(8)*ones(1,2),'Color','k','LineStyle',':')
                line(xl,yl(2)-subplot_start(8)+1.05/1.1*subplot_size(8)*ones(1,2),'Color','k','LineStyle',':')
                lw=1;
            else
                if subplot_posi(8)==subplot_posi(6)
                    maxmeanratemaval=maxrateval;
                else
                    maxmeanratemaval=maxmeanrateval;
                end
                lw=2.5;
            end
            plot(pratios_support,yl(2)-subplot_start(8)+(0.05+pmean_rates_ma/maxmeanratemaval)/1.1*subplot_size(8),['-',colors(4)],'LineWidth',lw)
        end

        if subplot_posi(9)>0                                                                          % A (pairwise average)
            maxbival=1;
            bilab=[0 0.5 1];
            yt=[yt yl(2)-subplot_start(9)+(0.05+bilab/maxbival)/1.1*subplot_size(9)];
            ytl=[ytl bilab];
            text(xl(1)-0.08*(xl(2)-xl(1)),yl(2)-subplot_start(9)+0.75/1.1*subplot_size(9),'A','Color','k','FontSize',11,'FontWeight','bold')
            line(xl,yl(2)-subplot_start(9)+0.05/1.1*subplot_size(9)*ones(1,2),'Color','k','LineStyle',':')
            line(xl,yl(2)-subplot_start(9)+1.05/1.1*subplot_size(9)*ones(1,2),'Color','k','LineStyle',':')
            if subplot_posi(10)~=subplot_posi(9)
                lw=2;
            else
                lw=1;
            end
            plot(pratios_support,yl(2)-subplot_start(9)+(0.05+ppa_ints/maxbival)/1.1*subplot_size(9),['-',colors(3)],'LineWidth',lw)

            if window_mode==2
                win_all_spikes2=all_spikes2(first_winspike_ind:last_winspike_ind);
                if wmin<all_spikes2(first_winspike_ind);
                    win_all_spikes2=[wmin win_all_spikes2];
                end
                if wmax>all_spikes2(last_winspike_ind);
                    win_all_spikes2=[win_all_spikes2 wmax];
                end
                win_pratios_support=sort([win_all_spikes2 win_all_spikes2(2:end-1)-dt]);
                win_ppa_ints=reshape([win_pa_ints; win_pa_ints],1,length(win_pa_ints)*2);
                plot(win_pratios_support,yl(2)-subplot_start(9)+(0.05+win_ppa_ints)/maxbival/1.1*subplot_size(9),'.-','Color',colors(5),'LineWidth',1)
                if patch_mode==1
                    patch([win_pratios_support(1) win_pratios_support win_pratios_support(end)],yl(2)-subplot_start(9)+(0.05+[0 win_ppa_ints 0])/1.1*subplot_size(9),colors(5))
                end
                line([wmin wmax],yl(2)-subplot_start(9)+(0.05+win_multi_pa_isi_dist/maxbival)/1.1*subplot_size(9)*ones(1,2),'LineWidth',1,'LineStyle','--','Color',colors(5))
            end
            % line([itmin itmax],yl(2)-subplot_start(9)+(0.05+0.5/maxbival)*subplot_size(9)/1.1*ones(1,2),'Color','k','LineStyle',':','LineWidth',2)
            line([itmin itmax],yl(2)-subplot_start(9)+(0.05+all_multi_pa_isi_dist/maxbival)/1.1*subplot_size(9)*ones(1,2),'Color',colors(3),'LineStyle','--')
        end

        if subplot_posi(10)>0                                                                          % A-MA (pairwise average)
            if subplot_posi(10)~=subplot_posi(9)
                maxbival=1;
                bilab=[0 0.5 1];
                yt=[yt yl(2)-subplot_start(10)+(0.05+bilab/maxbival)/1.1*subplot_size(10)];
                ytl=[ytl bilab];
                text(xl(1)-0.08*(xl(2)-xl(1)),yl(2)-subplot_start(10)+0.75/1.1*subplot_size(10),'A^*','Color','k','FontSize',11,'FontWeight','bold')
                line(xl,yl(2)-subplot_start(10)+0.05/1.1*subplot_size(10)*ones(1,2),'Color','k','LineStyle',':')
                line(xl,yl(2)-subplot_start(10)+1.05/1.1*subplot_size(10)*ones(1,2),'Color','k','LineStyle',':')
                line([itmin itmax],yl(2)-subplot_start(10)+(0.05+multi_pa_isi_dist/maxbival)/1.1*subplot_size(10)*ones(1,2),'Color','k','LineStyle','--')
                lw=1;
            else
                lw=2;
            end
            plot(pratios_support,yl(2)-subplot_start(10)+(0.05*subplot_size(10)+ppa_ints_ma/maxbival)/1.1*subplot_size(10),['-',colors(4)],'LineWidth',lw)
        end

        if subplot_posi(11)>0                                                                          % CV
            maxcvval=max([max(cv_ints(indy)) 1]);
            cvlab=0:maxcvval;
            yt=[yt yl(2)-subplot_start(11)+(0.05+cvlab/maxcvval)/1.1*subplot_size(11)];
            ytl=[ytl cvlab];
            text(xl(1)-0.08*(xl(2)-xl(1)),yl(2)-subplot_start(11)+0.75/1.1*subplot_size(11),'C_V','Color','k','FontSize',11,'FontWeight','bold')
            line(xl,yl(2)-subplot_start(11)+0.05/1.1*subplot_size(11)*ones(1,2),'Color','k','LineStyle',':')
            line(xl,yl(2)-subplot_start(11)+1.05/1.1*subplot_size(11)*ones(1,2),'Color','k','LineStyle',':')
            if subplot_posi(12)~=subplot_posi(11)
                lw=2;
            else
                lw=1;
            end
            plot(pratios_support,yl(2)-subplot_start(11)+(0.05+pcv_ints/maxcvval)/1.1*subplot_size(11),['-',colors(3)],'LineWidth',lw)

            if window_mode==2
                win_all_spikes2=all_spikes2(first_winspike_ind:last_winspike_ind);
                if wmin<all_spikes2(first_winspike_ind);
                    win_all_spikes2=[wmin win_all_spikes2];
                end
                if wmax>all_spikes2(last_winspike_ind);
                    win_all_spikes2=[win_all_spikes2 wmax];
                end
                win_pratios_support=sort([win_all_spikes2 win_all_spikes2(2:end-1)-dt]);
                win_pcv_ints=reshape([win_cv_ints; win_cv_ints],1,length(win_cv_ints)*2);
                plot(win_pratios_support,yl(2)-subplot_start(11)+(0.05+win_pcv_ints/maxcvval)/1.1*subplot_size(11),'.-','Color',colors(5),'LineWidth',1)
                if patch_mode==1
                    patch([win_pratios_support(1) win_pratios_support win_pratios_support(end)],yl(2)-subplot_start(11)+(0.05+[0 win_pcv_ints 0]/maxcvval)/1.1*subplot_size(11),colors(5))
                end
                line([wmin wmax],yl(2)-subplot_start(11)+(0.05+win_multi_cv_isi_dist/maxcvval)/1.1*subplot_size(11)*ones(1,2),'LineWidth',1,'LineStyle','--','Color',colors(5))
            end
            % line([itmin itmax],yl(2)-subplot_start(11)+(0.05+1/sqrt(2)/maxcvval)*subplot_size(11)/1.1*ones(1,2),'Color','k','LineStyle',':','LineWidth',2)
            line([itmin itmax],yl(2)-subplot_start(11)+(0.05+all_multi_cv_isi_dist/maxcvval)/1.1*subplot_size(11)*ones(1,2),'Color',colors(3),'LineStyle','--')
        end

        if subplot_posi(12)>0                                                                          % CV-MA
            if subplot_posi(12)~=subplot_posi(11)
                maxcvval=max([max(cv_ints(indy)) 1]);
                cvlab=0:maxcvval;
                yt=[yt yl(2)-subplot_start(12)+(0.05+cvlab/maxcvval)/1.1*subplot_size(12)];
                ytl=[ytl cvlab];
                line([itmin itmax],yl(2)-subplot_start(12)+(0.05+multi_cv_isi_dist/maxcvval)/1.1*subplot_size(12)*ones(1,2),'Color','k','LineStyle','--')
                % line([itmin itmax],yl(2)-subplot_start(12)+(0.05+1/maxcvval)*subplot_size(12)/1.1*ones(1,2),'Color','k','LineStyle',':','LineWidth',2)
                text(xl(1)-0.08*(xl(2)-xl(1)),yl(2)-subplot_start(12)+0.75/1.1*subplot_size(12),'C_V^*','Color','k','FontSize',11,'FontWeight','bold')
                line(xl,yl(2)-subplot_start(12)+0.05/1.1*subplot_size(12)*ones(1,2),'Color','k','LineStyle',':')
                line(xl,yl(2)-subplot_start(12)+1.05/1.1*subplot_size(12)*ones(1,2),'Color','k','LineStyle',':')
                lw=1;
            else
                lw=2;
            end
            plot(pratios_support,yl(2)-subplot_start(12)+(0.05+pcv_ints_ma/maxcvval)/1.1*subplot_size(12),['-',colors(4)],'LineWidth',lw)
        end

        if window_mode==2
            line(wmin*ones(1,2),yl,'Color',colors(5),'LineStyle','-.')
            line(wmax*ones(1,2),yl,'Color',colors(5),'LineStyle','--')
        end
        line(itmin*ones(1,2),yl,'Color','k','LineStyle',':')
        line(itmax*ones(1,2),yl,'Color','k','LineStyle',':')
        line(tmin*ones(1,2),yl,'Color','k','LineStyle','-.')
        line(tmax*ones(1,2),yl,'Color','k','LineStyle','-.')
        xlabel(['Time  ',timeunit_string],'Color','k','FontSize',12,'FontWeight','bold')
        title([title_string,'   ',comment_string,result_string],'Color','k','FontSize',14,'FontWeight','bold')
        [syt,syti]=sort(yt);
        sytl=ytl(syti);
        set(gca,'YTick',syt,'YTickLabel',sytl,'FontSize',10);
        set(gcf,'Color','w'); set(gca,'Color','w','Box','on');
        if print_mode                                                                    % Create postscript file
            set(gcf,'PaperOrientation','Landscape'); set(gcf,'PaperType', 'A4');
            set(gcf,'PaperUnits','Normalized','PaperPosition', [0 0 1.0 1.0]);
            psname=[filename,'_',num2str(plc),comment_string,'.ps'];
            print(gcf,'-dpsc',psname);
        end
    end
end
% uspikes(uspikes==987654321.123456789)=0;

%
% *************************************************************************
%
function z=f_moving_average(x,o)    % moving average past+future, inwards averaging at the edges (asymmetric)
if size(x,length(size(x)))==1
    x=x';
end
n=length(x);
if n<2*o+1
    o=fix((n-1)/2);
end
if o~=0
    for k=1:o
        z(k)=sum(x(k-(k-1):k+o))/((k-1)+1+o);
    end
    for k=1:2*o+1
        y(k,1:n-2*o)=x(k-1+(1:n-2*o));
    end
    z(o+1:n-o)=mean(y);
    for k=n-o+1:n
        z(k)=sum(x(k-o:k+(n-k)))/(o+1+(n-k));
    end
else
    z=x;
end
%
% *************************************************************************
%
function z=f_weighted_moving_average(y,x,o)    % past+future, inwards averaging at the edges (asymmetric), time weighted

if size(y,length(size(y)))==1
    y=y';
end
if size(x,length(size(x)))==1
    x=x';
end

n=length(y);
if n<2*o+1
    o=fix((n-1)/2);
end
if o~=0
    for k=1:o
        z(k)=sum(y(k-(k-1):k+o).*x(k-(k-1):k+o))/sum(x(k-(k-1):k+o));
    end
    for k=1:2*o+1
        y2(k,1:n-2*o)=y(k-1+(1:n-2*o));
        x2(k,1:n-2*o)=x(k-1+(1:n-2*o));
    end
    z(o+1:n-o)=sum(y2.*x2)./sum(x2,1);
    for k=n-o+1:n
        z(k)=sum(y(k-o:k+(n-k)).*x(k-o:k+(n-k)))/sum(x(k-o:k+(n-k)));
    end
else
    z=y;
end
%
% *************************************************************************
%
function lab=f_lab(vect,num_values,origin,cut)

if origin==1
    vect2=[0 vect];
else
    vect2=vect;
end
range=max(vect2)-min(vect2);
interval=range/num_values;
min_unit=fix(interval*10^(ceil(-log10(interval))))/10^(ceil(-log10(interval)));

xsep=round((range+min_unit)/num_values/min_unit)*min_unit;
lab=round(min(vect2)/min_unit)*min_unit+(xsep:xsep:(range+min_unit));

if cut==1
    lab=lab(lab>=min(vect2) & lab<=max(vect2));
end
%
% % *************************************************************************
%
function dt=f_get_dt(spikes)

vspikes=reshape(spikes,1,numel(spikes));

logs=-10:10;
c=0;
for lc=logs
    c=c+1;
    decades(c)=10^-lc;
    b=vspikes*decades(c);
    yes(c)=all(b==round(b));
    if yes(c)==0
        break;
    end
end

if c==1 || all(yes)
    dt=10^-8;
else
    dt=10^logs(c-1);
end