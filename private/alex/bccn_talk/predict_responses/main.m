%% load NSEM
mvpath='/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/';
my_movie=zeros(160,160,3600);
cnt=1;
for i=1:30
    load([mvpath, 'movie_chunk_', int2str(i), '.mat']);
    movie=movie(81:240,:,:);
%     movie=imresize(movie,2,'nearest');
    my_movie(:,:,cnt:cnt+119)=movie;
    cnt=cnt+120;
end
clear mvpath movie
my_movie = uint8(my_movie);
save('/Users/alexth/Desktop/BCCN2015/predictions/data/raw_movie', 'my_movie', '-v7.3');

%% prepare data

load('/Users/alexth/Desktop/BCCN2015/predictions/data/raw_movie')
my_movie=double(my_movie)/255;
my_movie = my_movie-mean(my_movie(:));

date = '2015-03-09-2';
date = '2015-08-17-1';

for ndf=4:-1:2
    [starun, nsemrun, wnrun, inputs, refresh, inputs_wn, refresh_wn, inputs_nsem] = load_stuff(date, ndf, my_movie);
    if size(inputs,2)>72000 % recording was long!!
        inputs = inputs(:,1:72000);
    end
    save(['/Users/alexth/Desktop/BCCN2015/predictions/data/data_ndf',int2str(ndf),'_', date],...
        'starun', 'nsemrun', 'wnrun', 'inputs',  'refresh', 'inputs_wn', 'refresh_wn', 'inputs_nsem', '-v7.3')
end


%% prepare asr
clear

trig_threshold = 0.84;
n_repeats = 20;

date = '2015-03-09-2';
date = '2015-08-17-1';

for ndf = 0:4
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/data_ndf',int2str(ndf),'_', date])
    
    clear asr asr_nsem sta sta_spikes wn_spikes nsem_spikes wn_trigs nsem_trigs asr_full inter_corr_wn inter_corr_nsem
    
    one_repeat_end = find(diff(wnrun.triggers)>trig_threshold,1);
    one_rep_wn = floor(wnrun.triggers(one_repeat_end)*1000/refresh); % duration of 1 repeat in frames
    one_repeat_end = find(diff(nsemrun.triggers)>trig_threshold,1);
    one_rep_nsem = floor(nsemrun.triggers(one_repeat_end)*1000/refresh); % duration of 1 repeat in frames
    
    
    cellIDs = get_cell_indices(starun, {1,2,3,4});
    good_cells = [];
    cnt = 1;
    for cellID = cellIDs
        % actual firing rate for WN and NSEM repeats
        
        [asr_full{cnt}, ~] = get_mean_spike_rate(wnrun, cellID, trig_threshold, n_repeats, 1000/8.3274);
        [asr{cnt}, inter_corr_wn(cnt)] = get_mean_spike_rate(wnrun, cellID, trig_threshold, n_repeats, 1000/refresh);
        [asr_nsem{cnt}, inter_corr_nsem(cnt)] = get_mean_spike_rate(nsemrun, cellID, trig_threshold, n_repeats, 1000/8.3274);
       
        sta_tmp = starun.stas.polarities{cellID}*squeeze(starun.stas.stas{cellID});
        sta_tmp = reshape(sta_tmp, [],30);
        sta_tmp = starun.stas.polarities{cellID}*sta_tmp(:,11:end);
        [~,max_loc] = find(sta_tmp==max(sta_tmp(:)));
        a = robust_std(sta_tmp(:,max_loc(1)));
        sta{cnt} = sta_tmp;
        
        sta_spikes{cnt} = starun.spikes{cellID}-starun.triggers(1);
        wn_spikes{cnt} = wnrun.spikes{cellID};
        nsem_spikes{cnt} = nsemrun.spikes{cellID};
        wn_trigs{cnt} = wnrun.triggers;
        nsem_trigs{cnt} = nsemrun.triggers;
        
        if ~isempty(find(sta_tmp(:,max_loc(1))>a*5, 1)) && numel(sta_spikes{cnt})>700 && ...
            length(find(sta_tmp(:,max_loc(1))>a*3.5))>3 && inter_corr_wn(cnt)>0.6 && inter_corr_nsem(cnt)>0.4
            good_cells = [good_cells cnt];
        end
        cnt = cnt+1;
    end
    figure
    plot(inter_corr_nsem)
    hold on
    plot(inter_corr_wn)
    drawnow
    
    save(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf',int2str(ndf),'_', date], ...
        'asr', 'asr_nsem', 'sta', 'sta_spikes', 'wn_spikes', 'nsem_spikes', 'wn_trigs', ...
        'nsem_trigs', 'good_cells', 'one_rep_wn', 'one_rep_nsem', 'asr_full', 'inter_corr_nsem', 'inter_corr_wn');
end


figure
std(inter_corr_nsem)


%% do modeling
clear

date = '2015-03-09-2';
date = '2015-08-17-1';

sta_params.length = 15;
sta_params.offset = 0;
fraction = 0.9;
trig_threshold = 0.84;
n_repeats = 20;
nbins = 100;

cnt = 1;
for ndf=2:4
    ndf
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf',int2str(ndf),'_', date]);
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/data_ndf',int2str(ndf),'_', date], 'refresh', 'inputs', 'inputs_wn', 'inputs_nsem')
    
    sta_rate = round(refresh/8.3);
    
    for i = good_cells
        
        %%%% estimate STA and nonlinearity from STA runs
        [~,max_loc] = find(sta{i}==max(sta{i}(:)));
        a = robust_std(sta{i}(:,max_loc(1)));
        my_stixels = find(sta{i}(:,max_loc(1))>a*3.5);
        
        train_spikes=ceil(sta_spikes{i}*1000/refresh)-one_rep_wn; % spikes in intervals
        train_spikes(train_spikes<1) = [];
        
        train_inputs = inputs(my_stixels,one_rep_wn+1:end); % in intervals
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA_coarse(train_inputs, train_spikes, fraction, sta_params);
        
        %%%%%%%%% predict response to WN
        
        gen_sig = conv2(inputs_wn(my_stixels,1:one_rep_wn), unbiased_sta,'valid');
        predicted_rate_wn{i, cnt} = convert_gs(gen_sig, gensig_bins, nonlinearity);
        
        r2_wn(i, cnt) = compute_r2(predicted_rate_wn{i, cnt}, asr{i}(sta_params.length:end));
        corr_wn(i, cnt) = compute_corr(predicted_rate_wn{i, cnt}, asr{i}(sta_params.length:end));
        
        %%%%%%%%% predict response to NSEM
        
        test_inputs_nsem = inputs_nsem(my_stixels,:); % assuming refresh every frame
        
        full_sta = interp1(1:sta_rate:sta_rate*sta_params.length, unbiased_sta', 1:sta_rate*sta_params.length, 'PCHIP')';
        full_sta = full_sta(:,1:end-sta_rate); % one sample point less
        
        gen_sig = conv2(test_inputs_nsem-mean(test_inputs_nsem(:)), full_sta,'valid');
        asr_tmp = asr_nsem{i}(size(full_sta,2):end);
        asr_tmp = asr_tmp(1:length(gen_sig));
        fit_params{i,cnt} = fit_nonlinearity(asr_tmp, gen_sig, nbins);
        
        predicted_rate_nsem{i, cnt} = fit_params{i,cnt}.scale*normcdf(gen_sig,fit_params{i,cnt}.mu, fit_params{i,cnt}.sigma);
        r2_nsem(i, cnt) = compute_r2(predicted_rate_nsem{i, cnt}, asr_tmp);
        corr_nsem(i, cnt) = compute_corr(predicted_rate_nsem{i, cnt}, asr_tmp);
   
    end
    cnt = cnt+1;
    
end


save(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_', date], ...
    'predicted_rate_wn', 'predicted_rate_nsem', 'r2_nsem', 'corr_nsem',...
    'r2_wn', 'corr_wn', 'fit_params')



%% NSEM: predict by other light level

date = '2015-08-17-1';

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf2_2015-08-17-1'], 'asr_nsem', 'good_cells');
ndf2 = asr_nsem;
gc2 = good_cells;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf3_2015-08-17-1'], 'asr_nsem', 'good_cells');
ndf3 = asr_nsem;
gc3 = good_cells;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf4_2015-08-17-1'], 'asr_nsem', 'good_cells');
ndf4 = asr_nsem;
gc4 = good_cells;

very_good_cells = intersect(intersect(gc2, gc3), gc4);
for i = very_good_cells
    %     r2_2(i, cnt) = compute_r2(ndf3, ndf2);
    corr_23(i) = compute_corr(ndf2{i}, ndf3{i});
    corr_24(i) = compute_corr(ndf2{i}, ndf4{i});
    corr_34(i) = compute_corr(ndf3{i}, ndf4{i});
end

save(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date], ...
    'corr_23', 'corr_24', 'corr_34','very_good_cells' )

clear

date = '2015-03-09-2';
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf0_2015-03-09-2'], 'asr_nsem', 'good_cells');
ndf0 = asr_nsem;
gc0 = good_cells;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf1_2015-03-09-2'], 'asr_nsem', 'good_cells');
ndf1 = asr_nsem;
gc1 = good_cells;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf2_2015-03-09-2'], 'asr_nsem', 'good_cells');
ndf2 = asr_nsem;
gc2 = good_cells;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf3_2015-03-09-2'], 'asr_nsem', 'good_cells');
ndf3 = asr_nsem;
gc3 = good_cells;
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf4_2015-03-09-2'], 'asr_nsem', 'good_cells');
ndf4 = asr_nsem;
gc4 = good_cells;

very_good_cells = intersect(intersect(gc2, gc3), gc1);
for i = very_good_cells
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
    'corr_04', 'corr_13', 'corr_03', 'corr_02','very_good_cells' )

%% plot NSEM inter-ndf correlations

clear

date = '2015-08-17-1';
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date])
figure
data = [corr_23(very_good_cells); corr_24(very_good_cells); corr_34(very_good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(very_good_cells)), 'x')
axis([0 4 0 1])
set(gca, 'xticklabel', {'2 and 3', '2 and 4', '3 and 4'})
xlabel('NDF')


clear
date = '2015-03-09-2';
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date]);
figure
% ndf4
subplot(2,3,1)
data = [corr_04(very_good_cells); corr_14(very_good_cells); corr_24(very_good_cells); corr_34(very_good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(very_good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'0', '1', '2', '3'})
title('NDF4')

% ndf3
subplot(2,3,2)
data = [corr_03(very_good_cells); corr_13(very_good_cells); corr_23(very_good_cells); corr_34(very_good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(very_good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'0', '1', '2', '4'})
title('NDF3')

% ndf2
subplot(2,3,3)
data = [corr_02(very_good_cells); corr_12(very_good_cells); corr_23(very_good_cells); corr_24(very_good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(very_good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'0', '1', '3', '4'})
title('NDF2')

% ndf1
subplot(2,3,4)
data = [corr_01(very_good_cells); corr_12(very_good_cells); corr_13(very_good_cells); corr_14(very_good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(very_good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'0', '2', '3', '4'})
title('NDF1')

% ndf0
subplot(2,3,5)
data = [corr_01(very_good_cells); corr_02(very_good_cells); corr_03(very_good_cells); corr_04(very_good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(very_good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'1', '2', '3', '4'})
title('NDF0')

% all by previous
subplot(2,3,6)
data = [corr_01(very_good_cells); corr_12(very_good_cells); corr_23(very_good_cells); corr_34(very_good_cells)]';
bar(mean(data))
hold on
errorbar(mean(data), std(data)/sqrt(numel(very_good_cells)), 'x')
axis([0 5 0 1])
set(gca, 'xticklabel', {'01', '12', '23', '34'})



%% plot correlations and r2: modeling
clear
date = '2015-03-09-2';

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_', date])
figure
for i=1:5
    subplot(2,3,i)
    tmp = corr_wn(:,i);
    tmp1 = corr_nsem(:,i);
    
    plot(tmp(tmp~=0), tmp1(tmp1~=0), 'x')
    
    tt(1,i) = mean(tmp(tmp~=0));
    tt(2,i) = std(tmp(tmp~=0))/sqrt(nnz(tmp));
    tt1(1,i) = mean(tmp1(tmp1~=0));
    tt1(2,i) = std(tmp1(tmp1~=0))/sqrt(nnz(tmp1));
    title(int2str(i-1))
    axis([0 1 0 1])
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

% load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date]);
% figure
% for i=1:5
%     subplot(2,3,i)
%     clear tmp
%     for j=very_good_cells
%         tmp(1,j) = fit_params{j,i}.scale;
%         tmp(2,j) = fit_params{j,i}.mu;
%         tmp(3,j) = fit_params{j,i}.sigma;
%     end
%     plot(tmp(3,:))
% end


clear
date = '2015-08-17-1';

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_', date])
figure
for i=1:3
    subplot(1,3,i)
    tmp = corr_wn(:,i);
    tmp1 = corr_nsem(:,i);
    
    plot(tmp(tmp~=0), tmp1(tmp1~=0), 'x')
    
    tt(1,i) = mean(tmp(tmp~=0));
    tt(2,i) = std(tmp(tmp~=0))/sqrt(nnz(tmp));
    tt1(1,i) = mean(tmp1(tmp1~=0));
    tt1(2,i) = std(tmp1(tmp1~=0))/sqrt(nnz(tmp1));
    title(int2str(i-1))
    axis([0 1 0 1])
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

%% model vs inter-ndf predictions
clear
date = '2015-03-09-2';

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_', date])
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date]);

figure
% ndf 4 and 3
subplot(2,3,1)
model_ndf = corr_nsem(very_good_cells,5);
inter_ndf = corr_34(very_good_cells);
plot(inter_ndf, model_ndf, 'x')
title('NDF 4')
axis([0 1 0 1])
xlabel('NDF 3, data')
ylabel('NDF 4, model')

% ndf 3 and 2
subplot(2,3,2)
model_ndf = corr_nsem(very_good_cells,4);
inter_ndf = corr_24(very_good_cells);
plot(inter_ndf, model_ndf, 'x')
title('NDF 3')
axis([0 1 0 1])
xlabel('NDF 4, data')
ylabel('NDF 3, model')

% ndf 2 and 1
subplot(2,3,3)
model_ndf = corr_nsem(very_good_cells,3);
inter_ndf = corr_23(very_good_cells);
plot(inter_ndf, model_ndf, 'x')
title('NDF 2')
axis([0 1 0 1])
xlabel('NDF 3, data')
ylabel('NDF 2, model')

% ndf 1 and 0
subplot(2,3,4)
model_ndf = corr_nsem(very_good_cells,2);
inter_ndf = corr_01(very_good_cells);
plot(inter_ndf, model_ndf, 'x')
title('NDF 1')
axis([0 1 0 1])
xlabel('NDF 0, data')
ylabel('NDF 1, model')

% ndf 0 and 1
subplot(2,3,5)
model_ndf = corr_nsem(very_good_cells,1);
inter_ndf = corr_01(very_good_cells);
plot(inter_ndf, model_ndf, 'x')
title('NDF 0')
axis([0 1 0 1])
xlabel('NDF 1, data')
ylabel('NDF 0, model')



date = '2015-08-17-1';

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_', date])
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date]);

figure
% ndf 4 and 3
subplot(1,3,1)
model_ndf = corr_nsem(very_good_cells,3);
inter_ndf = corr_34(very_good_cells);
plot(inter_ndf, model_ndf, 'x')
title('NDF 4')
axis([0 1 0 1])
xlabel('NDF 3, data')
ylabel('NDF 4, model')

% ndf 3 and 2
subplot(1,3,2)
model_ndf = corr_nsem(very_good_cells,2);
inter_ndf = corr_23(very_good_cells);
plot(inter_ndf, model_ndf, 'x')
title('NDF 3')
axis([0 1 0 1])
xlabel('NDF 2, data')
ylabel('NDF 3, model')

% ndf 2 and 3
subplot(1,3,3)
model_ndf = corr_nsem(very_good_cells,1);
inter_ndf = corr_23(very_good_cells);
plot(inter_ndf, model_ndf, 'x')
title('NDF 2')
axis([0 1 0 1])
xlabel('NDF 3, data')
ylabel('NDF 2, model')


%% plot rasters: 2015-03-09-2
clear
date = '2015-03-09-2';
close all
ndf1 = 2;
ndf2 = 4;

load(['/Users/alexth/Desktop/BCCN2015/predictions/data/predictions_', date])
load(['/Users/alexth/Desktop/BCCN2015/predictions/data/nsem_correlations_', date], 'very_good_cells');


for cellID = very_good_cells(2:10);
    
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf',int2str(ndf1),'_', date]);
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/data_ndf',int2str(ndf1),'_', date], 'refresh')
    sta_rate = round(refresh/8.3);
    one_frame = refresh/sta_rate;
    otst = 14*sta_rate;

    dt = 1;
    fr = predicted_rate_nsem{cellID,ndf1+1};
    nTrials = 14;
    nBins = length(fr);
    
    figure
    set(gcf, 'position', [-1887         461        1842         619]);
    clear h
    
    % prediction light level 1
    h(1) = subplot(4,1,1);
    spikeMat = rand(nTrials, nBins) < repmat(fr*dt,14,1);
    hold off
    for j=1:14
        tmp = find(spikeMat(j,:));
        plot(tmp+otst,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
        hold on
    end
    title(['Model, ',date,', cellID ',int2str(cellID),', NDF', int2str(ndf1), ...
        ', corr = ', num2str(corr_nsem(cellID,ndf1+1)), ', r2 = ', num2str(r2_nsem(cellID,ndf1+1))]);
    
    % actual data light level 1
    h(2) = subplot(4,1,2);
    spikes=nsem_spikes{cellID};
    trigs = nsem_trigs{cellID};
    myTrigs=[0 find(diff(trigs)>0.84)'];
    splitSpikes=cell(14,1);
    for j=6:19
        tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
            - trigs(myTrigs(j)+1);
        splitSpikes{j-5}=ceil(tmp*1000/one_frame);
    end
    hold off
    tt = 0;
    for j=1:14
        plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
        tt = tt + numel(splitSpikes{j});
        hold on
    end
    title(['Data, ',date,', cellID ',int2str(cellID),', NDF', int2str(ndf1)]);
    
    
    
    
    fr = predicted_rate_nsem{cellID,ndf2+1};
    nBins = length(fr);
    
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/asr_ndf',int2str(ndf2),'_', date]);
    load(['/Users/alexth/Desktop/BCCN2015/predictions/data/data_ndf',int2str(ndf2),'_', date], 'refresh')
    sta_rate = round(refresh/8.3);
    one_frame = refresh/sta_rate;
    otst = 14*sta_rate;
    
    % actual data light level 2
    h(3) = subplot(4,1,3);
    spikes=nsem_spikes{cellID};
    trigs = nsem_trigs{cellID};
    myTrigs=[0 find(diff(trigs)>0.84)'];
    splitSpikes=cell(14,1);
    for j=6:19
        tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
            - trigs(myTrigs(j)+1);
        splitSpikes{j-5}=ceil(tmp*1000/one_frame);
    end
    hold off
    tt = 0;
    for j=1:14
        plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
        tt = tt + numel(splitSpikes{j});
        hold on
    end
    title(['Data, ',date,', cellID ',int2str(cellID),', NDF', int2str(ndf2)]);
    
    
    % prediction light level 2
    h(4) = subplot(4,1,4);
    spikeMat = rand(nTrials, nBins) < repmat(fr*dt,14,1);
    hold off
    for j=1:14
        tmp = find(spikeMat(j,:));
        plot(tmp+otst,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
        hold on
    end
    title(['Model, ',date,', cellID ',int2str(cellID),', NDF', int2str(ndf2), ...
        ', corr = ', num2str(corr_nsem(cellID,ndf2+1)), ', r2 = ', num2str(r2_nsem(cellID,ndf2+1))]);
    
end


linkaxes(h,'xy')
axis([100 3600 0 0.8])


%% plot rasters: example

cellIDs = get_cell_indices(starun{1}, {1});
date = '2015-08-17-1';
ndfs = '34';

plot_nsem = 1;
plot_wn = 0;

for cellID = cellIDs
    
    clear r2_wn_wrongND corr_wn_wrongND r2_nsem_wrongND corr_nsem_wrongND
    
    
    
    asr{1,cnt} = get_mean_spike_rate(wnrun{1}, cellID, trig_threshold, n_repeats, 1000/refresh_rep(1));
    asr{2,cnt} = get_mean_spike_rate(wnrun{2}, cellID, trig_threshold, n_repeats, 1000/refresh_rep(2));
    asr_nsem{1,cnt} = get_mean_spike_rate(nsemrun{1}, cellID, trig_threshold, n_repeats, 1000/refresh_rep(1));
    asr_nsem{2,cnt} = get_mean_spike_rate(nsemrun{2}, cellID, trig_threshold, n_repeats, 1000/refresh_rep(2));
    
    r2_wn_wrongND(1) = compute_r2(asr{2,cnt}, asr{1,cnt});
    r2_wn_wrongND(2) = compute_r2(asr{1,cnt}, asr{2,cnt});
    corr_wn_wrongND(1) = compute_corr(asr{2,cnt}, asr{1,cnt});
    r2_nsem_wrongND(1) = compute_r2(asr_nsem{2,cnt}, asr_nsem{1,cnt});
    r2_nsem_wrongND(2) = compute_r2(asr_nsem{1,cnt}, asr_nsem{2,cnt});
    corr_nsem_wrongND(1) = compute_corr(asr_nsem{2,cnt}, asr_nsem{1,cnt});
    
    
    % light level 1
    predicted_rate = cell(1,2);
    predicted_rate_nsem = cell(1,2);
    corr_wn = zeros(1,2);
    r2_wn = zeros(1,2);
    r2_nsem = zeros(1,2);
    corr_nsem = zeros(1,2);
    
    for run_number = 1:2
        
        sta = squeeze(starun{run_number}.stas.stas{cellID});
        sta = reshape(sta, [],30);
        sta = starun{run_number}.stas.polarities{cellID}*sta(:,11:end);
        [~,max_loc] = find(sta==max(sta(:)));
        a = robust_std(sta(:,max_loc(1)))*4;
        my_stixels = find(sta(:,max_loc(1))>a);
        
        
        train_inputs = inputs{run_number}(my_stixels,601:end); % in intervals
        train_spikes=ceil((starun{run_number}.spikes{cellID}-starun{run_number}.triggers(1))*1000/refresh(run_number))-600; % spikes in intervals
        train_spikes(train_spikes<1) = [];
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA_coarse(train_inputs, train_spikes, fraction, sta_params);
        
        % white noise
        test_inputs_wn = inputs_rep{run_number}(my_stixels,1:one_rep_wn);
        gen_sig = conv2(test_inputs_wn, unbiased_sta,'valid');
        predicted_rate{run_number} = convert_gs(gen_sig, gensig_bins, nonlinearity);
        
        corr_wn(run_number) = compute_corr(predicted_rate{run_number},asr{run_number}(sta_params.length:end));
        r2_wn(run_number) = compute_r2(predicted_rate{run_number},asr{run_number}(sta_params.length:end));
        
        % nsem
        test_inputs_nsem = inputs_nsem{run_number}(my_stixels,1:one_rep_nsem); % assuming refresh every frame
        gen_sig = conv2(test_inputs_nsem, unbiased_sta,'valid');
        asr_tmp = asr_nsem{run_number}(sta_params.length:end);
        asr_tmp = asr_tmp(1:length(gen_sig));
        [new_nl, gensig_bins] = get_nonlinearity(asr_tmp, gen_sig, nbins);
        predicted_rate_nsem{run_number} = convert_gs(gen_sig, gensig_bins, new_nl);
        
        r2_nsem(run_number) = compute_r2(predicted_rate_nsem{run_number},asr_tmp);
        corr_nsem(run_number) = compute_corr(predicted_rate_nsem{run_number},asr_tmp);
    end
    
    
    % plot white noise
    
    if plot_wn
        figure
        set(gcf, 'position', [-1887         461        1842         619]);
        clear h
        
        % prediction light level 1
        h(1) = subplot(4,1,1);
        fr = predicted_rate{1}; % spikes per frames
        dt = 0.83; % frame
        nBins = length(fr); % 10 ms spike train
        nTrials = 14; % number of simulations
        spikeMat = rand(nTrials, nBins) < repmat(fr*dt,14,1);
        a = sum(spikeMat(:));
        hold off
        for j=1:14
            tmp = find(spikeMat(j,:));
            plot(tmp+14,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
            hold on
        end
        title(['Model, ',date,', cellID ',int2str(cellID),', NDF', ndfs(1), ...
            ', corr = ', num2str(corr_wn(1)), ', r2 = ', num2str(r2_wn(1))]);
        
        % actual data light level 1
        h(2) = subplot(4,1,2);
        spikes=wnrun{1}.spikes{cellID};
        trigs = wnrun{1}.triggers;
        myTrigs=[0 find(diff(trigs)>0.84)'];
        splitSpikes=cell(14,1);
        for j=6:19
            tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
                - trigs(myTrigs(j)+1);
            splitSpikes{j-5}=ceil(tmp*1000/refresh_rep(1));
        end
        hold on
        tt = 0;
        for j=1:14
            plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
            tt = tt + numel(splitSpikes{j});
        end
        title(['Data, ',date,', cellID ',int2str(cellID),', NDF', ndfs(1), ...
            ', corr = ', num2str(corr_wn_wrongND(1)), ', r2 = ', num2str(r2_wn_wrongND(1))]);
        
        
        % actual data light level 2
        h(3) = subplot(4,1,3);
        spikes=wnrun{2}.spikes{cellID};
        trigs = wnrun{2}.triggers;
        myTrigs=[0 find(diff(trigs)>0.84)'];
        splitSpikes=cell(14,1);
        for j=6:19
            tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
                - trigs(myTrigs(j)+1);
            splitSpikes{j-5}=ceil(tmp*1000/refresh_rep(1));
        end
        hold on
        tt = 0;
        for j=1:14
            plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
            tt = tt + numel(splitSpikes{j});
        end
        title(['Data, ',date,', cellID ',int2str(cellID),', NDF', ndfs(2), ...
            ', corr = ', num2str(corr_wn_wrongND(1)), ', r2 = ', num2str(r2_wn_wrongND(2))]);
        
        % prediction light level 2
        h(4) = subplot(4,1,4);
        fr = predicted_rate{2}; % spikes per frames
        dt = 0.83; % frame
        nBins = length(fr);
        nTrials = 14; % number of simulations
        spikeMat = rand(nTrials, nBins) < repmat(fr*dt,14,1);
        a = sum(spikeMat(:));
        hold off
        for j=1:14
            tmp = find(spikeMat(j,:));
            plot(tmp+14,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
            hold on
        end
        title(['Model, ',date,', cellID ',int2str(cellID),', NDF', ndfs(2), ...
            ', corr = ', num2str(corr_wn(2)), ', r2 = ', num2str(r2_wn(2))]);
        
        
        linkaxes(h,'xy')
        
        axis([20 length(predicted_rate{2}) 0 0.8])
    end
    
    % plot nsem
    if plot_nsem
        figure
        set(gcf, 'position', [-1887         461        1842         619]);
        clear h
        
        % prediction light level 1
        h(1) = subplot(4,1,1);
        fr = predicted_rate_nsem{1}; % spikes per frames
        dt = 0.83; % frame
        nBins = length(fr); % 10 ms spike train
        nTrials = 14; % number of simulations
        spikeMat = rand(nTrials, nBins) < repmat(fr*dt,14,1);
        a = sum(spikeMat(:));
        hold off
        for j=1:14
            tmp = find(spikeMat(j,:));
            plot(tmp+14,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
            hold on
        end
        title(['Model, ',date,', cellID ',int2str(cellID),', NDF', ndfs(1), ...
            ', corr = ', num2str(corr_nsem(1)), ', r2 = ', num2str(r2_nsem(1))]);
        
        % actual data light level 1
        h(2) = subplot(4,1,2);
        spikes=nsemrun{1}.spikes{cellID};
        trigs = nsemrun{1}.triggers;
        myTrigs=[0 find(diff(trigs)>0.84)'];
        splitSpikes=cell(14,1);
        for j=6:19
            tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
                - trigs(myTrigs(j)+1);
            splitSpikes{j-5}=ceil(tmp*1000/refresh_rep(1));
        end
        hold on
        tt = 0;
        for j=1:14
            plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
            tt = tt + numel(splitSpikes{j});
        end
        title(['Data, ',date,', cellID ',int2str(cellID),', NDF', ndfs(1), ...
            ', corr = ', num2str(corr_nsem_wrongND(1)), ', r2 = ', num2str(r2_nsem_wrongND(1))]);
        
        
        % actual data light level 2
        h(3) = subplot(4,1,3);
        spikes=nsemrun{2}.spikes{cellID};
        trigs = nsemrun{2}.triggers;
        myTrigs=[0 find(diff(trigs)>0.84)'];
        splitSpikes=cell(14,1);
        for j=6:19
            tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
                - trigs(myTrigs(j)+1);
            splitSpikes{j-5}=ceil(tmp*1000/refresh_rep(1));
        end
        hold on
        tt = 0;
        for j=1:14
            plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
            tt = tt + numel(splitSpikes{j});
        end
        title(['Data, ',date,', cellID ',int2str(cellID),', NDF', ndfs(2), ...
            ', corr = ', num2str(corr_nsem_wrongND(1)), ', r2 = ', num2str(r2_nsem_wrongND(2))]);
        
        % prediction light level 2
        h(4) = subplot(4,1,4);
        fr = predicted_rate_nsem{2}; % spikes per frames
        dt = 0.83; % frame
        nBins = length(fr);
        nTrials = 14; % number of simulations
        spikeMat = rand(nTrials, nBins) < repmat(fr*dt,14,1);
        a = sum(spikeMat(:));
        hold off
        for j=1:14
            tmp = find(spikeMat(j,:));
            plot(tmp+14,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
            hold on
        end
        title(['Model, ',date,', cellID ',int2str(cellID),', NDF', ndfs(2), ...
            ', corr = ', num2str(corr_nsem(2)), ', r2 = ', num2str(r2_nsem(2))]);
        
        linkaxes(h,'xy')
        axis([20 length(predicted_rate{2}) 0 0.8])
    end
    
end
