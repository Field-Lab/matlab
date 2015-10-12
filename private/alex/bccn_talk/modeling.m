%% prepare data
clear

load('/Users/alexth/Desktop/BCCN2015/predictions/data/raw_movie')
my_movie=double(my_movie)/255;
my_movie = my_movie-mean(my_movie(:));

date = '2015-03-09-2';
date = '2015-08-17-1';

for ndf=4:-1:2
    [inputs, inputs_wn, inputs_nsem, meta, spikes, sta] = load_stuff(date, ndf, my_movie);

    save(['/Users/alexth/Desktop/BCCN2015/predictions/data/inputs_ndf',int2str(ndf),'_', date],...
        'inputs', 'inputs_wn', 'inputs_nsem', '-v7.3');
 
    save(['/Users/alexth/Desktop/BCCN2015/predictions/data/spikes_ndf',int2str(ndf),'_', date],...
        'meta', 'spikes', 'sta');
end

%% prepare asr
clear

trig_threshold = 0.84;
n_repeats = 20;
fps = 1000/8.3274;

date = '2015-03-09-2';
date = '2015-08-17-1';

meth = 'count';

for ndf = 0:4
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/spikes_ndf',int2str(ndf),'_', date])
    
    clear asr asr_nsem
    for i=1:size(sta,2)
        [asr_nsem{i}, inter_corr_nsem(i)] = get_mean_spike_rate(spikes.nsem{i}, meta.triggers.nsem, trig_threshold, n_repeats, fps, meth);
        [asr{i}, inter_corr(i)] = get_mean_spike_rate(spikes.wn{i}, meta.triggers.wn, trig_threshold, n_repeats, fps, meth);
    end
    
    save(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf',int2str(ndf),'_', date], ...
        'asr', 'asr_nsem', 'inter_corr_nsem', 'inter_corr')
end
    
%% estimate properties: STA and nonlinearity
clear

date = '2015-03-09-2';
date = '2015-08-17-1';

trig_threshold = 0.84;
sta_params.length = 15;
sta_params.offset = 0;
fraction = 0.5;

for ndf = 2:4
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/spikes_ndf',int2str(ndf),'_', date])
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/inputs_ndf',int2str(ndf),'_', date], 'inputs')
    
    one_repeat_end = find(diff(meta.triggers.wn)>trig_threshold,1);
    one_rep_wn = floor(meta.triggers.wn(one_repeat_end)*1000/meta.refresh); % duration of 1 repeat in frames
    one_repeat_end = find(diff(meta.triggers.nsem)>trig_threshold,1);
    one_rep_nsem = floor(meta.triggers.nsem(one_repeat_end)*1000/meta.refresh); % duration of 1 repeat in frames
   
    clear unbiased_sta fit_params my_stixels gof
    for i=1:size(sta,2)
        i        
        tmp = sta(:, i);
        a = robust_std(tmp(:));
        if a>0
            tmp_stixels = find(abs(tmp)>a*2);
            if ~isempty(tmp_stixels)
                train_spikes = ceil(spikes.sta{i}*1000/meta.refresh-one_rep_wn);
                train_spikes(train_spikes<1) = [];
                train_inputs = inputs(tmp_stixels,one_rep_wn+1:end)*0.96-0.48;
                
                [~, ~, ~, stix_sign]=unbiased_STA_coarse(train_inputs, train_spikes, fraction, sta_params);
                
                my_stixels{i} = tmp_stixels(stix_sign>fraction);
                if ~isempty(my_stixels{i})
                    train_inputs = inputs(my_stixels{i},one_rep_wn+1:end)*0.96-0.48;
                    [unbiased_sta{i}, gensig_bins, nonlinearity, stix_sign]=unbiased_STA_coarse(train_inputs, train_spikes, fraction, sta_params);
                    gen_sig = conv2(train_inputs, unbiased_sta{i},'valid');
                    train_spikes(train_spikes>size(train_inputs,2)) = [];
                    train_spikes(train_spikes<sta_params.length) = [];
                    spikes_tmp = train_spikes;
                    spike_rate=zeros(size(train_inputs,2),1);
                    while ~isempty(spikes_tmp)
                        [~, ia, ~] = unique(spikes_tmp);
                        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
                        spikes_tmp(ia)=[];
                    end
                    fitres = fit(gensig_bins(1:end-1)',nonlinearity, 'scale*normcdf(x,mu,sigma)','StartPoint',[1, max(nonlinearity), 1.5], 'Lower', [0.1 0.1 0.1], 'Upper', [15 15 15]);
                    [fitres1, gof{i}] = fit(gen_sig',spike_rate(sta_params.length:end), 'scale*normcdf(x,mu,sigma)', 'StartPoint', [fitres.mu, fitres.scale, fitres.sigma], 'Lower', [0.1 0.1 0.1], 'Upper', [15 15 15]);
                    fit_params{i}.scale = fitres1.scale;
                    fit_params{i}.mu = fitres1.mu;
                    fit_params{i}.sigma = fitres1.sigma;
                end
            end
        end
    end
    
    save(['/Users/alexth/Desktop/BCCN2015/predictions/data/sta_ndf',int2str(ndf),'_', date], ...
        'unbiased_sta', 'fit_params', 'my_stixels', 'gof')
    
end



%% do modeling
clear

date = '2015-03-09-2';
date = '2015-08-17-1';

sta_params.length = 15;
sta_params.offset = 0;
trig_threshold = 0.84;
nbins = 100;

for ndf=2:4
    ndf
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/sta_ndf',int2str(ndf),'_', date]);
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/inputs_ndf',int2str(ndf),'_', date], 'inputs_wn', 'inputs_nsem')
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/spikes_ndf',int2str(ndf),'_', date], 'meta')
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf',int2str(ndf),'_', date], 'asr', 'asr_nsem')
    
    one_repeat_end = find(diff(meta.triggers.wn)>trig_threshold,1);
    one_rep_wn = floor(meta.triggers.wn(one_repeat_end)*1000/meta.refresh); % duration of 1 repeat in frames
    one_repeat_end = find(diff(meta.triggers.nsem)>trig_threshold,1);
    one_rep_nsem = floor(meta.triggers.nsem(one_repeat_end)*1000/meta.refresh); % duration of 1 repeat in frames
    
    sta_rate = round(meta.refresh/8.3);
    
    clear predicted_rate_wn predicted_rate_nsem actual_rate_wn actual_rate_nsem
    clear r2_wn corr_wn r2_nsem corr_nsem fit_params_nsem
    
    for i = 1:size(asr,2)
        
        if ~isempty(fit_params{i}) %&& fit_params{i}.sigma <2 && gof{i}.rsquare>0.1
            
            %%%%%%%%% predict response to WN
            inputs = inputs_wn(my_stixels{i},1:one_rep_wn)*0.96-0.48;
            gen_sig = conv2(inputs, unbiased_sta{i},'valid');
            predicted_rate_wn{i} = fit_params{i}.scale*normcdf(gen_sig,fit_params{i}.mu, fit_params{i}.sigma);
            actual_rate = sum(reshape(asr{i}(1:sta_rate*one_rep_wn),sta_rate, []));
            actual_rate_wn{i} = actual_rate(sta_params.length:end);
            
            r2_wn(i) = compute_r2(predicted_rate_wn{i}, actual_rate_wn{i});
            corr_wn(i) = compute_corr(predicted_rate_wn{i}, actual_rate_wn{i});
            
            
            %%%%%%%%% predict response to NSEM
            full_sta = interp1(1:sta_rate:sta_rate*sta_params.length, unbiased_sta{i}', 1:sta_rate*sta_params.length, 'PCHIP')';
            full_sta = full_sta(:,1:end-sta_rate); % one sample point less
            
            inputs = inputs_nsem(my_stixels{i},:); % assuming refresh every frame
            
            gen_sig = conv2(inputs, full_sta,'valid');
            asr_tmp = asr_nsem{i}(size(full_sta,2):end);
            actual_rate_nsem{i} = asr_tmp(1:length(gen_sig));
            fit_params_nsem{i} = fit_nonlinearity(actual_rate_nsem{i}, gen_sig, nbins);
            
            predicted_rate_nsem{i} = fit_params_nsem{i}.scale*normcdf(gen_sig,fit_params_nsem{i}.mu, fit_params_nsem{i}.sigma);
            r2_nsem(i) = compute_r2(predicted_rate_nsem{i}, actual_rate_nsem{i});
            corr_nsem(i) = compute_corr(predicted_rate_nsem{i}, actual_rate_nsem{i});
            
        end
    end
    save(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf),'_', date], ...
        'predicted_rate_wn', 'predicted_rate_nsem', 'r2_nsem', 'corr_nsem',...
        'r2_wn', 'corr_wn', 'fit_params_nsem', 'actual_rate_wn', 'actual_rate_nsem')
    
end

%% get cell types

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data011/data011';
datarun = load_data(path2data);
datarun = load_params(datarun,'verbose',1);
for i=1:length(datarun.cell_ids)
    [type_name{i},type_index(i)]=find_cell_type(datarun, datarun.cell_ids(i));
end

save(['/Users/alexth/Desktop/BCCN2015/predictions/data/cell_types_', date], ...
    'type_name', 'type_index')

%% NSEM: predict by other light level
clear

date = '2015-08-17-1';

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf2_2015-08-17-1'], 'asr_nsem', 'inter_corr_nsem');
ndf2 = asr_nsem;
inter_corr_nsem2 = inter_corr_nsem;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf3_2015-08-17-1'], 'asr_nsem', 'inter_corr_nsem');
ndf3 = asr_nsem;
inter_corr_nsem3 = inter_corr_nsem;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf4_2015-08-17-1'], 'asr_nsem', 'inter_corr_nsem');
ndf4 = asr_nsem;
inter_corr_nsem4 = inter_corr_nsem;

for i = 1:size(ndf2,2)
    corr_23(i) = compute_corr(ndf2{i}, ndf3{i});
    corr_24(i) = compute_corr(ndf2{i}, ndf4{i});
    corr_34(i) = compute_corr(ndf3{i}, ndf4{i});
end

save(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date], ...
    'corr_23', 'corr_24', 'corr_34', 'inter_corr_nsem2', 'inter_corr_nsem3', 'inter_corr_nsem4')

clear

date = '2015-03-09-2';
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf0_2015-03-09-2'],'asr_nsem', 'inter_corr_nsem');
ndf0 = asr_nsem;
inter_corr_nsem0 = inter_corr_nsem;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf1_2015-03-09-2'],'asr_nsem', 'inter_corr_nsem');
ndf1 = asr_nsem;
inter_corr_nsem1 = inter_corr_nsem;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf2_2015-03-09-2'],'asr_nsem', 'inter_corr_nsem');
ndf2 = asr_nsem;
inter_corr_nsem2 = inter_corr_nsem;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf3_2015-03-09-2'],'asr_nsem', 'inter_corr_nsem');
ndf3 = asr_nsem;
inter_corr_nsem3 = inter_corr_nsem;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf4_2015-03-09-2'],'asr_nsem', 'inter_corr_nsem');
ndf4 = asr_nsem;
inter_corr_nsem4 = inter_corr_nsem;

for i = 1:size(ndf2,2)
    corr_01(i) = compute_corr(ndf0{i}, ndf1{i});
    corr_12(i) = compute_corr(ndf1{i}, ndf2{i});
    corr_23(i) = compute_corr(ndf2{i}, ndf3{i});
    corr_34(i) = compute_corr(ndf3{i}, ndf4{i});
    
    corr_24(i) = compute_corr(ndf2{i}, ndf4{i});
    corr_14(i) = compute_corr(ndf1{i}, ndf4{i});
    corr_04(i) = compute_corr(ndf0{i}, ndf4{i});
    
    corr_13(i) = compute_corr(ndf1{i}, ndf3{i});
    corr_03(i) = compute_corr(ndf0{i}, ndf3{i});
    
    corr_02(i) = compute_corr(ndf0{i}, ndf2{i});
end

save(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date], ...
    'corr_01', 'corr_12', 'corr_23', 'corr_34', 'corr_24', 'corr_14',...
    'corr_04', 'corr_13', 'corr_03', 'corr_02','inter_corr_nsem0', 'inter_corr_nsem1', ...
    'inter_corr_nsem2', 'inter_corr_nsem3', 'inter_corr_nsem4')

%% plot NSEM inter-ndf correlations

clear

date = '2015-08-17-1';
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date])

figure
plot(inter_corr_nsem2)
hold on
plot(inter_corr_nsem3)
plot(inter_corr_nsem4)

good_cells = find(inter_corr_nsem2>0.5 & inter_corr_nsem3>0.5 & inter_corr_nsem4>0.5);


figure
data = [corr_23(good_cells); corr_24(good_cells); corr_34(good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(size(data,1)), 'x')
axis([0 4 0 1])
set(gca, 'xticklabel', {'2 and 3', '2 and 4', '3 and 4'})
xlabel('NDF')


clear
date = '2015-03-09-2';
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date]);

figure
plot(inter_corr_nsem0)
hold on
plot(inter_corr_nsem1)
plot(inter_corr_nsem2)
plot(inter_corr_nsem3)
plot(inter_corr_nsem4)

good_cells = find(inter_corr_nsem0>0.5 & inter_corr_nsem1>0.5 & inter_corr_nsem2>0.5 & inter_corr_nsem3>0.5 & inter_corr_nsem4>0.5);

figure
% ndf4
subplot(2,3,1)
data = [corr_04(good_cells); corr_14(good_cells); corr_24(good_cells); corr_34(good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'0', '1', '2', '3'})
title('NDF4')

% ndf3
subplot(2,3,2)
data = [corr_03(good_cells); corr_13(good_cells); corr_23(good_cells); corr_34(good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'0', '1', '2', '4'})
title('NDF3')

% ndf2
subplot(2,3,3)
data = [corr_02(good_cells); corr_12(good_cells); corr_23(good_cells); corr_24(good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'0', '1', '3', '4'})
title('NDF2')

% ndf1
subplot(2,3,4)
data = [corr_01(good_cells); corr_12(good_cells); corr_13(good_cells); corr_14(good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'0', '2', '3', '4'})
title('NDF1')

% ndf0
subplot(2,3,5)
data = [corr_01(good_cells); corr_02(good_cells); corr_03(good_cells); corr_04(good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'1', '2', '3', '4'})
title('NDF0')

% all by previous
subplot(2,3,6)
data = [corr_01(good_cells); corr_12(good_cells); corr_23(good_cells); corr_34(good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'01', '12', '23', '34'})



%% plot correlations and r2: modeling
clear
date = '2015-03-09-2';

figure
for ndf=0:4
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf),'_', date])
    subplot(2,3,ndf+1)
    
    tmp = corr_wn(corr_wn~=0);
    tmp1 = corr_nsem(corr_nsem~=0);
    plot(tmp, tmp1, 'x')
    
    tt(1, ndf+1) = mean(tmp);
    tt(2, ndf+1) = std(tmp)/sqrt(numel(tmp));
    tt1(1, ndf+1) = mean(tmp1);
    tt1(2, ndf+1) = std(tmp1)/sqrt(numel(tmp1));
    title(int2str(ndf))
    axis([0 1 0 1])
    line([0 1], [0 1], 'color', 'k')
    xlabel('white noise')
    ylabel('nsem')
end

figure
subplot(1,2,1)
bar(tt(1,:))
hold on
errorbar(tt(1,:), tt(2,:), 'x')
set(gca, 'xticklabel', {'0', '1', '2', '3', '4'})
title('white noise')
axis([0 6 0 1])
subplot(1,2,2)
bar(tt1(1,:))
hold on
errorbar(tt1(1,:), tt1(2,:), 'x')
set(gca, 'xticklabel', {'0', '1', '2', '3', '4'})
title('nsem')
axis([0 6 0 1])


clear
date = '2015-08-17-1';

figure
for ndf=2:4
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf),'_', date])
    subplot(1,3,ndf-1)
    
    tmp = corr_wn(corr_wn~=0);
    tmp1 = corr_nsem(corr_nsem~=0);
    plot(tmp, tmp1, 'x')
    
    tt(1, ndf-1) = mean(tmp);
    tt(2, ndf-1) = std(tmp)/sqrt(numel(tmp));
    tt1(1, ndf-1) = mean(tmp1);
    tt1(2, ndf-1) = std(tmp1)/sqrt(numel(tmp1));
    title(int2str(ndf))
    axis([0 1 0 1])
    line([0 1], [0 1], 'color', 'k')
    xlabel('white noise')
    ylabel('nsem')
end

figure
subplot(1,2,1)
bar(tt(1,:))
hold on
errorbar(tt(1,:), tt(2,:), 'x')
set(gca, 'xticklabel', { '2', '3', '4'})
title('white noise')
axis([0 4 0 1])
subplot(1,2,2)
bar(tt1(1,:))
hold on
errorbar(tt1(1,:), tt1(2,:), 'x')
set(gca, 'xticklabel', {'2', '3', '4'})
title('nsem')
axis([0 4 0 1])

%% model vs inter-ndf predictions: NSEM
clear
date = '2015-03-09-2';

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date]);
ndf2 = [1 0 3 4 3];

figure
for ndf = 0:4
%     load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf2(ndf+1)),'_', date])    
%     second_corr = max(corr_nsem, 0.5);
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf),'_', date])    
%     primary_corr = max(corr_nsem, 0.5);
    
    data_corr = eval(['corr_', int2str(min(ndf, ndf2(ndf+1))), int2str(max(ndf, ndf2(ndf+1)))]);
    good_cells = eval(['find(inter_corr_nsem',int2str(ndf),'>0.5 & inter_corr_nsem',int2str(ndf2(ndf+1)),'>0.5);']);
    good_cells = intersect(good_cells, find(corr_nsem>0));
%     good_cells(good_cells>length(corr_nsem))=[];
    subplot(2,3,ndf+1)
    hold on
    for k=1:length(unique(type_index))
        if sum(type_index==k)>0 && ~strcmp(type_name(find(type_index==k,1)), 'duplicates')
            my_cells = intersect(good_cells, find(type_index==k));
            model_ndf = corr_nsem(my_cells);
            inter_ndf = data_corr(my_cells);
            plot(inter_ndf, model_ndf, 'x')
        end
    end
    
    
%     model_ndf = corr_nsem(good_cells);
%     inter_ndf = data_corr(good_cells);
%     plot(inter_ndf, model_ndf, 'x')
    title(['NDF ', int2str(ndf), ', n = ', int2str(length(good_cells))])
    axis([0 1 0 1])
    xlabel(['NDF ', int2str(ndf2(ndf+1)),', DATA'])
    ylabel(['NDF ', int2str(ndf),', MODEL'])
    line([0, 1], [0, 1], 'color', 'k')
end


% summary plot
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/cell_types_', date])
clear my_corrs my_corrs_std
for ndf = 0:4
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf),'_', date]) 
    good_cells = eval(['find(inter_corr_nsem',int2str(ndf),'>0.5);']);
    for k=1:4
%         my_cells = intersect(good_cells, find(type_index==k));
        my_cells = find(type_index==k);
        my_corrs(k, ndf+1) = mean(corr_nsem(my_cells));
        my_corrs_std(k, ndf+1) = std(corr_nsem(my_cells))/sqrt(length(my_cells));
    end
end
figure
for i=1:4
    subplot(2,2,i)
    bar(my_corrs(i,:));
    axis([0 6 0 1])
end
figure
for i=1:5
    subplot(1,5,i)    
    bar(my_corrs(:,i));
    hold on
    errorbar(my_corrs(:,i), my_corrs_std(:,i), 'x')
    axis([0 5 0 1])
end


% nonlinearities across cell types and NDFs
date = '2015-03-09-2';
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/cell_types_', date])
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date]);
x = -5:0.1:5;
figure
for ndf = 0:4
%         load(['/Users/alexth/Desktop/BCCN2015/predictions/data/sta_ndf',int2str(ndf),'_', date],'fit_params')
        load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf),'_', date], 'fit_params_nsem') 
    good_cells = eval(['find(inter_corr_nsem',int2str(ndf),'>0.5);']);
    
    for k=1:4
        figure
        hold on
        my_cells = intersect(good_cells, find(type_index==k));
        a = [];
%         my_cells = find(type_index==k);
        for i=my_cells
            if ~isempty(fit_params{i})
%                 plot(x,fit_params_nsem{i}.scale*normcdf(x,fit_params_nsem{i}.mu, fit_params_nsem{i}.sigma));
                a = [a;fit_params{i}.sigma ];
%                 a = [a; fit_params{i}.scale*normcdf(x,fit_params{i}.mu, fit_params{i}.sigma)];
            end
        end
        plot(a)
    end
end





load(['/Users/alexth/Desktop/BCCN2015/predictions/data/cell_types_', date])
figure
my_type = cell(1);
for ndf = 0:4
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf),'_', date])    
    prim_inter = eval(['inter_corr_nsem',int2str(ndf)']);
    subplot(2,3,ndf+1)
    hold on
    for k=1:length(unique(type_index))
        if sum(type_index==k)>0 && ~strcmp(type_name(find(type_index==k,1)), 'duplicates')
            plot(prim_inter(type_index==k), corr_nsem(type_index==k), 'x')
            if ndf==0
                my_type = [my_type type_name(find(type_index==k,1))];
            end
        end        
    end
    if ndf==4        
        legend(my_type(2:end), 'location', 'best')
    end
    title(['NDF ', int2str(ndf)])
    axis([0 1 0 1])
    xlabel('Inter-data correlation')
    ylabel('Model correlation')
    line([0, 1], [0, 1], 'color', 'k')
end





load(['/Users/alexth/Desktop/BCCN2015/predictions/data/cell_types_', date])
figure
my_type = cell(1);
for ndf = 0:4
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf),'_', date])    
    prim_inter = eval(['inter_corr_nsem',int2str(ndf)']);
    
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf',int2str(ndf),'_',date],'asr_nsem');
    t = [];
    for j = 1:212
        t(j) = mean(asr_nsem{j});
    end
%     t = t/max(t)*max(corr_nsem);
    
    subplot(2,3,ndf+1)
    hold on
    for k=1:length(unique(type_index))
        if sum(type_index==k)>0 && ~strcmp(type_name(find(type_index==k,1)), 'duplicates')
             plot(t(type_index==k), prim_inter(type_index==k), 'x')
%             plot(prim_inter(type_index==k), corr_nsem(type_index==k), 'x')
            if ndf==0
                my_type = [my_type type_name(find(type_index==k,1))];
            end
        end        
    end
    if ndf==4        
        legend(my_type(2:end), 'location', 'best')
    end
    title(['NDF ', int2str(ndf)])
%     axis([0 1 0 1])
    xlabel('std in data')
    ylabel('prim_inter correlation')
%     line([0, 1], [0, 1], 'color', 'k')
end




date = '2015-08-17-1';

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date]);
figure
ndf2 = [3 4 3];
for ndf = 2:4
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf),'_', date])    
    
    data_corr = eval(['corr_', int2str(min(ndf, ndf2(ndf-1))), int2str(max(ndf, ndf2(ndf-1)))]);
    good_cells = eval(['find(inter_corr_nsem',int2str(ndf),'>0.5 & inter_corr_nsem',int2str(ndf2(ndf-1)),'>0.5);']);
    good_cells(good_cells>length(corr_nsem))=[];
    subplot(2,3,ndf-1)
    model_ndf = corr_nsem(good_cells);
    inter_ndf = data_corr(good_cells);
    plot(inter_ndf, model_ndf, 'x')
    title(['NDF ', int2str(ndf)])
    axis([0 1 0 1])
    xlabel(['NDF ', int2str(ndf2(ndf-1)),', DATA'])
    ylabel(['NDF ', int2str(ndf),', MODEL'])
    line([0, 1], [0, 1], 'color', 'k')
end

%% model vs inter-ndf predictions: WN
clear
date = '2015-03-09-2';

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf1_',date],'asr');
ndf1 = asr;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf2_',date],'asr');
ndf2 = asr;

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf1_', date])  
corr_wn1 = corr_wn;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf2_', date])  
corr_wn2 = corr_wn;

for i=1:212
    inter_corr_wn(i) = corr(ndf1{i}, ndf2{i});
end
figure
plot(inter_corr_wn, corr_wn2, '*')
line([0, 1], [0, 1], 'color', 'k')


clear
date = '2015-08-17-1';

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf3_',date],'asr');
ndf3 = asr;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf4_',date],'asr');
ndf4 = asr;

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf3_', date])  
corr_wn3 = corr_wn;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf4_', date])  
corr_wn4 = corr_wn;

for i=1:146
    inter_corr_wn(i) = corr(ndf3{i}, ndf4{i});
end
figure
plot(inter_corr_wn, corr_wn4, '*')
line([0, 1], [0, 1], 'color', 'k')




%% plot rasters: NSEM
clear
date = '2015-03-09-2';
date = '2015-08-17-1';
close all
ndf1 = 1;
ndf2 = 4;

dt = 1;
nTrials = 14;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf1),'_', date], 'predicted_rate_nsem', 'corr_nsem')   
psr1 = predicted_rate_nsem;
corr_nsem1 = corr_nsem;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/spikes_ndf',int2str(ndf1),'_', date], 'spikes', 'meta')
spikes1 = spikes.nsem;
trigs1 = meta.triggers.nsem;
refresh1 = meta.refresh;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf2),'_', date], 'predicted_rate_nsem', 'corr_nsem')   
psr2 = predicted_rate_nsem;
corr_nsem2 = corr_nsem;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/spikes_ndf',int2str(ndf2),'_', date], 'spikes', 'meta')
spikes2 = spikes.nsem;
trigs2 = meta.triggers.nsem;
refresh2 = meta.refresh;

sta_rate1 = round(refresh1/8.3);
otst1 = 14*sta_rate1;
sta_rate2 = round(refresh2/8.3);
otst2 = 14*sta_rate2;


one_frame = refresh1/sta_rate1;

for i = 69
    
    
    if ~isempty(psr1{i}) && ~isempty(psr2{i})
        
        figure
        set(gcf, 'position', [-1887         461        1842         619]);
        clear h
        
        % prediction light level 1
        h(1) = subplot(4,1,1);
        fr = psr1{i};
        nBins = length(fr);
        spikeMat = rand(nTrials, nBins) < repmat(fr*dt,14,1);
        hold off
        for j=1:14
            tmp = find(spikeMat(j,:));
            plot(tmp+otst1,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.r', 'markersize',0.1)
            hold on
        end
        title(['Model, ',date,', cell ',int2str(i),', NDF', int2str(ndf1), ...
            ', corr = ', num2str(corr_nsem1(i))]);
        
        % actual data light level 1
        h(2) = subplot(4,1,2);
        spikes=spikes1{i};
        myTrigs=[0 find(diff(trigs1)>0.84)'];
        hold off
        for j=6:19
            tmp=spikes(spikes>=trigs1(myTrigs(j)+1) & spikes<trigs1(myTrigs(j+1)))...
                - trigs1(myTrigs(j)+1);
            tmp=ceil(tmp*1000/one_frame);
            plot(tmp,zeros(length(tmp),1)+0.8/13*(j-5)-0.8/13, '.k', 'markersize',0.1)
            hold on
        end
        title(['Data, NDF', int2str(ndf1)]);
        
        
        % actual data light level 2
        h(3) = subplot(4,1,3);
        spikes=spikes2{i};
        myTrigs=[0 find(diff(trigs2)>0.84)'];
        hold off
        for j=6:19
            tmp=spikes(spikes>=trigs2(myTrigs(j)+1) & spikes<trigs2(myTrigs(j+1)))...
                - trigs2(myTrigs(j)+1);
            tmp=ceil(tmp*1000/one_frame);
            plot(tmp,zeros(length(tmp),1)+0.8/13*(j-5)-0.8/13, '.k', 'markersize',0.1)
            hold on
        end
        title(['Data, NDF', int2str(ndf2)]);
        
        
        % prediction light level 2
        h(4) = subplot(4,1,4);
        fr = psr2{i};
        nBins = length(fr);
        spikeMat = rand(nTrials, nBins) < repmat(fr*dt,14,1);
        hold off
        for j=1:14
            tmp = find(spikeMat(j,:));
            plot(tmp+otst2,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.r', 'markersize',0.1)
            hold on
        end
        title(['Model, ',date,', cell ',int2str(i),', NDF', int2str(ndf2), ...
            ', corr = ', num2str(corr_nsem2(i))]);
        
        
        linkaxes(h,'xy')
        axis([150 3600 0 0.8])
    end
    

end

for i=1:4
    subplot(h(i))
    axis off
end

%% plot rasters: white noise
clear
date = '2015-03-09-2';
date = '2015-08-17-1';
close all
ndf1 = 1;
ndf2 = 2;

dt = 1;
nTrials = 14;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf1),'_', date], 'predicted_rate_wn', 'corr_wn')   
psr1 = predicted_rate_wn;
corr_nsem1 = corr_wn;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/spikes_ndf',int2str(ndf1),'_', date], 'spikes', 'meta')
spikes1 = spikes.wn;
trigs1 = meta.triggers.wn;
refresh1 = meta.refresh;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf2),'_', date], 'predicted_rate_wn', 'corr_wn')   
psr2 = predicted_rate_wn;
corr_nsem2 = corr_wn;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/spikes_ndf',int2str(ndf2),'_', date], 'spikes', 'meta')
spikes2 = spikes.wn;
trigs2 = meta.triggers.wn;
refresh2 = meta.refresh;

sta_rate1 = round(refresh1/8.3);
otst1 = 14;
sta_rate2 = round(refresh2/8.3);
otst2 = 14;


one_frame = refresh1/sta_rate1;

for i = 99
    
    
    if ~isempty(psr1{i}) && ~isempty(psr2{i})
        
        figure
        set(gcf, 'position', [-1887         461        1842         619]);
        clear h
        
        % prediction light level 1
        h(1) = subplot(4,1,1);
        fr = psr1{i};
        nBins = length(fr);
        spikeMat = rand(nTrials, nBins) < repmat(fr*dt,14,1);
        hold off
        for j=1:14
            tmp = find(spikeMat(j,:))*sta_rate1;
            tmp = tmp+(rand(1, size(tmp,2)))*sta_rate1/2;
            plot(tmp+otst1*sta_rate1,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.r', 'markersize',0.1)
            hold on
        end
        title(['Model, ',date,', cell ',int2str(i),', NDF', int2str(ndf1), ...
            ', corr = ', num2str(corr_nsem1(i))]);
        
        % actual data light level 1
        h(2) = subplot(4,1,2);
        spikes=spikes1{i};
        myTrigs=[0 find(diff(trigs1)>0.84)'];
        hold off
        for j=6:19
            tmp=spikes(spikes>=trigs1(myTrigs(j)+1) & spikes<trigs1(myTrigs(j+1)))...
                - trigs1(myTrigs(j)+1);
            tmp=ceil(tmp*1000/one_frame);
            plot(tmp,zeros(length(tmp),1)+0.8/13*(j-5)-0.8/13, '.k', 'markersize',0.1)
            hold on
        end
        title(['Data, NDF', int2str(ndf1)]);
        
        
        % actual data light level 2
        h(3) = subplot(4,1,3);
        spikes=spikes2{i};
        myTrigs=[0 find(diff(trigs2)>0.84)'];
        hold off
        for j=6:19
            tmp=spikes(spikes>=trigs2(myTrigs(j)+1) & spikes<trigs2(myTrigs(j+1)))...
                - trigs2(myTrigs(j)+1);
            tmp=ceil(tmp*1000/one_frame);
            plot(tmp,zeros(length(tmp),1)+0.8/13*(j-5)-0.8/13, '.k', 'markersize',0.1)
            hold on
        end
        title(['Data, NDF', int2str(ndf2)]);
        
        % prediction light level 2
        h(4) = subplot(4,1,4);
        fr = psr2{i};
        nBins = length(fr);
        spikeMat = rand(nTrials, nBins) < repmat(fr*dt,14,1);
        hold off
        for j=1:14
            tmp = find(spikeMat(j,:))*sta_rate2;
            tmp = tmp+(rand(1, size(tmp,2))-0.5)*sta_rate2;
            plot(tmp+otst2*sta_rate2,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.r', 'markersize',0.1)
            hold on
        end
        title(['Model, ',date,', cell ',int2str(i),', NDF', int2str(ndf2), ...
            ', corr = ', num2str(corr_nsem2(i))]);
        
        linkaxes(h,'xy')
        axis([120 3400 0 0.8])
    end    

end

for i=1:4
    subplot(h(i))
    axis([150 1700 0 0.8])
end



%% RMS
clear
date = '2015-03-09-2';
% date = '2015-08-17-1';
close all

cag_data = zeros(4,5);
cag = zeros(4,5);
cag_data_std = zeros(4,5);
cag_std = zeros(4,5);
ddd =  zeros(4,5);
ddd_std = zeros(4,5);
kk=1;
for ndf = 0:4
    
    dt = 1;
    nTrials = 14;
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_ndf',int2str(ndf),'_', date], 'predicted_rate_nsem', 'corr_nsem')
    psr1 = predicted_rate_nsem;
    corr_nsem1 = corr_nsem;
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/spikes_ndf',int2str(ndf),'_', date], 'spikes', 'meta')
    spikes1 = spikes.nsem;
    trigs1 = meta.triggers.nsem;
    refresh1 = meta.refresh;
    
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/cell_types_', date])
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date]);
    
    otst1 = 14*round(refresh1/8.3);
    
    trig_threshold = 0.84;
    n_repeats = 20;
    fps = 1000/8.3274;
    
    sr=120;
    sig=40;
    st=10000/sr*6.2/(60*sig);
    time_list=-3.1:st:3.1;
    kern=zeros(1,length(time_list));
    for i=1:length(time_list)
        kern(i)=250/sig*exp((1-time_list(i)^2)/2);
    end
    figure
    plot(time_list, kern)
%     
%     figure
    
    for uu = 1:4
        myc = find(type_index == uu);
%         subplot(2,2,uu)
        cnt = 1;
        mean_corr = [];
        mean_corr_data = [];
        for i = myc%69%myc
            
            
            if ~isempty(psr1{i}) %&& eval(['inter_corr_nsem',int2str(ndf),'(i)>0.5'])
                
                % model firing rate
                fr = psr1{i};
                
                % data firing rates
                spikes=spikes1{i};
                repeats = get_repeats_spike_rate(spikes, trigs1, trig_threshold, n_repeats, fps);
                for k=1:20
                    repeats(k,:) = conv(repeats(k,:),kern,'same');
                end
                fr = conv(fr,kern,'same');
                t = [];
                for k=-15:15
                    t(k+16) = corr(mean(repeats(:,otst1-k:end-k-122))', fr');
                end
                [~, b] = max(t);
                lag = -15+b-1;
                
                my_repeats = repeats(:,otst1-lag:end-lag-122);
                
                % smooth!!!
                fr = fr(150:3400);
                my_repeats = my_repeats(:, 150:3400);
                corr(fr', mean(my_repeats)');

%                         subplot(2,1,1)
% figure
% hold on
% plot(mean(my_repeats(5:19,:)), 'b', 'linewidth',2)
% plot(fr, 'r', 'linewidth',2)
% title(int2str(ndf))
% compute_corr(mean(my_repeats(5:19,:)), fr)
                %
                %         mean(p)
                %         corr(fr',mean(my_repeats(1:k-1,:))')
                
%                 figure
%                 hold on
%                 tmp = my_repeats(5:19,:);
%                 my_trial = tmp(6,:);
%                 tmp(6,:) = [];
%                 plot(my_trial, 'r', 'linewidth',2)
%                 hold on
%                 plot(mean(tmp),'b', 'linewidth',2)
%                 compute_corr(my_trial, mean(tmp))
                
%                 figure
%                 hold on
%                 tmp = my_repeats(5:19,:);
%                 my_trial = tmp(6,:);
%                 tmp(6,:) = [];
%                 plot(my_trial, 'b', 'linewidth',2)
%                 hold on
%                 plot(fr,'k', 'linewidth',2)
%                 compute_corr(my_trial, fr)
                
                
                model_corr = [];data_corr=[];
                model_r2 = [];data_r2=[];
                model_rms = [];data_rms=[];
                for k=5:19
                    tmp = my_repeats(5:19,:);
                    tmp(k-4,:) = [];
                    model_corr(k-4) = compute_corr(fr,my_repeats(k,:));
                    data_corr(k-4) = compute_corr(my_repeats(k,:), mean(tmp));
                    model_r2(k-4) = compute_r2(fr,my_repeats(k,:));
                    data_r2(k-4) = compute_r2(my_repeats(k,:), mean(tmp));
                    model_rms(k-4) = rms(my_repeats(k,:)-fr);
                    data_rms(k-4) = rms(my_repeats(k,:)-mean(tmp));
                end
                
%                 figure
%                 plot(t, 'r*-')
%                 hold on
%                 plot(p, 'b*-')
%                 axis([0 16 0.5 1])
%                 mean(t)/mean(p)
%                 
                
%                 hold on
                if (nanmean(data_corr))>0.3 && (nanmean(data_r2))>0.3 && (nanmean(data_rms))<3.5
                    mean_corr(cnt) = nanmean(model_corr);
                    mean_corr_data(cnt) = nanmean(data_corr);
                    mean_r2(cnt) = nanmean(model_r2);
                    mean_r2_data(cnt) = nanmean(data_r2);
                    mean_rms(cnt) = nanmean(model_rms);
                    mean_rms_data(cnt) = nanmean(data_rms);
%                     plot(p,t,'x');
                else
                    mean_corr(cnt) = -100;
                    mean_corr_data(cnt) = -100;
                end
                
%                 if mean_corr(cnt)>0.85
%                     figure(kk)
%                     hold on
%                     plot(mean(my_repeats), 'k')
%                     plot(fr, 'r', 'linewidth',2)
%                     title([int2str(i), '  NDF', int2str(ndf)])
%                     kk = kk+1;
%                     drawnow
%                 end
%                 figure(40+ndf)
                cnt = cnt+1;
                %         subplot(2,1,2)
                %         hold on
                %         plot(p, 'k')
                %         plot(t, 'r')
                %         axis([0 15 0 1])
                %         title(['Cell ', int2str(i), ', ', type_name{i}, ', corr ', num2str(corr_nsem1(i))]);
                
            end
            
        end
%         xlabel('data')
%         ylabel('model')
%         line([0 1], [0 1], 'color', 'k')
%         title(type_name{myc(1)});
        gg = mean_corr_data>0.3;
        cag(uu, ndf+1) = mean(mean_corr(gg));
        cag_data(uu, ndf+1) = mean(mean_corr_data(gg));
        corr_trend_mean(uu, ndf+1) = mean(mean_corr(gg)./mean_corr_data(gg));
        corr_trend_std(uu, ndf+1) = std(mean_corr(gg)./mean_corr_data(gg))/sqrt(sum(gg));
        
        r2_trend_mean(uu, ndf+1) = mean(mean_r2(gg)./mean_r2_data(gg));
        r2_trend_std(uu, ndf+1) = std(mean_r2(gg)./mean_r2_data(gg))/sqrt(sum(gg));
        
        
        rms_trend_mean(uu, ndf+1) = mean(mean_rms_data(gg)./mean_rms(gg));
        rms_trend_std(uu, ndf+1) = std(mean_rms_data(gg)./mean_rms(gg))/sqrt(sum(gg));
        
       
        cag_std(uu, ndf+1) = std(mean_corr(gg))/sqrt(sum(gg));
        cag_data_std(uu, ndf+1) = std(mean_corr_data(gg))/sqrt(sum(gg));
    end
end

figure

subplot(1,3,1)
for i=1:4
    plot(corr_trend_mean(i, end:-1:1))
    hold on
end
legend('ONp', 'off p', 'onm', 'off p', 'location', 'best')
for i=1:4
errorbar(corr_trend_mean(i, end:-1:1), corr_trend_std(i, end:-1:1), '.')

end
axis([0 6 0 1])

subplot(1,3,2)
for i=1:4
    plot(r2_trend_mean(i, end:-1:1))
    hold on
end
legend('ONp', 'off p', 'onm', 'off p', 'location', 'best')
for i=1:4
errorbar(r2_trend_mean(i, end:-1:1), r2_trend_std(i, end:-1:1), '.')

end
axis([0 6 0 1])

subplot(1,3,3)
for i=1:4
    plot(rms_trend_mean(i, end:-1:1))
    hold on
end
legend('ONp', 'off p', 'onm', 'off p', 'location', 'best')
for i=1:4
errorbar(rms_trend_mean(i, end:-1:1), rms_trend_std(i, end:-1:1), '.')

end
axis([0 6 0 1])








figure
plot(mean_corr_data, mean_corr, 'o', 'markersize',10)
hold on
plot(mean_corr_data(12), mean_corr(12), 'o', 'markersize',10)
line([0 1], [0 1], 'color', 'r')


figure
for i=1:5
    subplot(2,5,i)
    bar(cag(:,i))
    hold on
    errorbar(cag(:,i), cag_std(:,i), '.r')
    axis([0 5 0 1])
    subplot(2,5,i+5)
    bar(cag_data(:,i))
    hold on
    errorbar(cag_data(:,i), cag_data_std(:,i), '.r')
    axis([0 5 0 1])
end



figure
for i=1:5
    subplot(1,5,i)
    bar(ddd(:,i))
    axis([0 5 0 1])
end


