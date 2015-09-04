function [mean_unbiased_sta, gensig_bins, mean_nonlinearity, stix_sign]=unbiased_STA_coarse(inputs, spikes, fraction, sta_params)
% this function takes certain fraction of cone inputs and calculates STA
% on it. Then calculates nonlinearity on remaining inputs. Repeats several
% times choosing different snippets of inputs for STA and NL. Takes average
% of STAs and average of NLs for each cone, resulting in an unbiased STA
% and NL. 
% Spikes have to start from time 0 as inputs do.
% fraction - what fraction of inputs to take for each sta (0.5, 0.9, etc)
% first, last portion of inputs will be taken for nonlinearity
% then, next adjacent portion


%parameters
bin_nonl = 100; % number of bins for nonlinearity

spikes(spikes>size(inputs,2)) = []; 
spikes(spikes<sta_params.length) = [];
% get spike rate
spikes_tmp = spikes;
spike_rate=zeros(size(inputs,2),1);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp
n_stixels = size(inputs,1);

% common STA to estimate bins for nonlinearity
spikes_tmp = spikes;
spikes_tmp(spikes_tmp<sta_params.length) = []; 
sta=zeros(n_stixels,sta_params.length);
nspikes = numel(spikes_tmp);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    for j=1:sta_params.length
        sta(:,sta_params.length-j+1) = sta(:,sta_params.length-j+1)...
            + sum(inputs(:,spikes_tmp(ia) - sta_params.length + j + sta_params.offset),2);
    end
    spikes_tmp(ia)=[];
end
sta = sta/nspikes;


% generator signal bins
    
tmp_spike_rate = spike_rate(sta_params.length:end);

gen_sig = conv2(inputs, sta,'valid');

% gen_sig=gen_sig/max(gen_sig);
bin_size = ceil(size(gen_sig,2)/bin_nonl);
tmp = sort(gen_sig);
gensig_bins = [tmp(1:bin_size:end) tmp(end)*1.001];

my_nl=zeros(bin_nonl,1);
% mean firing rate in each bin of generator signal
for k=1:bin_nonl
    my_nl(k)=mean(tmp_spike_rate(gen_sig>=gensig_bins(k) & gen_sig<gensig_bins(k+1)));
end
    

total_stas = round(1/(1 - fraction)); % number of partial STAs to be calculated
nonl_fraction =floor((1-fraction)*size(inputs,2)); % size of data to exclude from partial STA (for nonlinearity)

all_nl = zeros(bin_nonl, total_stas); % stores all partial nonlinearities
all_partial_stas = zeros(n_stixels,sta_params.length, total_stas);

% figure
for sta_cnt=1:total_stas
    nonlin_spikes = nonl_fraction*[sta_cnt-1 sta_cnt];
    
    
    spikes_for_sta = spikes;
    % remove spikes which wil be taken for nonlinearity
    spikes_for_sta(spikes_for_sta > nonlin_spikes(1) & spikes_for_sta < nonlin_spikes(2)) = [];
    nspikes = numel(spikes_for_sta);
    
    % partial inputs and spike rate
    partial_inputs = inputs(:,(nonlin_spikes(1)+1):nonlin_spikes(2));    
    partial_spike_rate = spike_rate((nonlin_spikes(1)+sta_params.length):nonlin_spikes(2));

    % partial sta
    partial_sta=zeros(n_stixels,sta_params.length);
    
    while ~isempty(spikes_for_sta)
        [~, ia, ~] = unique(spikes_for_sta);
        for j=1:sta_params.length
            partial_sta(:,sta_params.length-j+1) = partial_sta(:,sta_params.length-j+1)...
                + sum(inputs(:,spikes_for_sta(ia) - sta_params.length + j + sta_params.offset),2);
        end
        spikes_for_sta(ia)=[];
    end
    partial_sta = partial_sta/nspikes;
    all_partial_stas(:,:,sta_cnt) = partial_sta;
    
    gen_sig = conv2(partial_inputs, partial_sta,'valid');
%     gen_sig=gen_sig/max(gen_sig);
               
    tmp = gensig_bins;
    my_nl=zeros(bin_nonl,1);
    % mean firing rate in each bin of generator signal
    for k=1:bin_nonl
        my_nl(k)=mean(partial_spike_rate(gen_sig>=tmp(k) & gen_sig<tmp(k+1)));
    end
    all_nl(:,sta_cnt) = my_nl;
    
end
% set(ax, 'YLim',[min(all_nl(:)), max(all_nl(:))])

stix_sign = zeros(1,size(sta,1));
for i=1:size(sta,1)
    a = squeeze(all_partial_stas(i,:,:));
    tmp = triu(corr(a),1);
    tmp = tmp(:);
    tmp(tmp==0) = [];
    stix_sign(i) = mean(tmp);
end

mean_unbiased_sta = nanmean(all_partial_stas,3);
mean_nonlinearity = nanmean(all_nl,2);

if isnan(mean_nonlinearity(1)); mean_nonlinearity(1) = 0; end
for k=2:length(mean_nonlinearity)
    if isnan(mean_nonlinearity(k))
        mean_nonlinearity(k) = mean_nonlinearity(k-1);
    end
end

% figure
% for i=1:n_stixels
%     tmp = gensig_bins(i,:);
%     ax(i) = subplot(4,5,i);
%     plot(tmp(1:end-1),mean_nonlinearity(i,:))
% end
% set(ax, 'YLim',[min(mean_nonlinearity(:)), max(mean_nonlinearity(:))])



