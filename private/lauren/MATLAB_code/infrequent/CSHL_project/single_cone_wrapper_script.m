clear all
close all

datarun = [];
datapath = 'C:\Documents and Settings\Lauren Hruby\Desktop\plantain';

datarun = import_single_cone_data(datarun, datapath);
%note: distance unit ~= 4 microns

%scale all distances to ~= microns

datarun.cone_centers = 4*datarun.cone_centers;
datarun.rgc_COMs = 4*datarun.rgc_COMs;

%thresholding minimum cone weight for each cell
for ii = 1:length(datarun.cell_ids)
    thresh = 0.15*max(abs(datarun.cone_weights(:,ii)));
    datarun.cone_weights(datarun.cone_weights(:,ii)<thresh, ii) = 0;
end


includeNeuron = false(length(datarun.cell_ids), 1);
minConePos = [1000 1000];
conesIncluded = zeros(size(datarun.cone_centers,1),1); %binary vector specifying which cones are connected to chosen cells
for ii = 1:length(datarun.cell_ids)
    if datarun.cell_types(ii) == 3 || datarun.cell_types(ii) == 4;
        center = datarun.rgc_COMs(ii,:);
        if center(1) < 200*4 && center(1) > 140*4 && center(2) > 155*4 && center(2) < 200*4
            includeNeuron(ii) = true;
            minConePos(1) = min(minConePos(1), min(datarun.cone_centers(datarun.cone_weights(:,ii)>0,1)));
            minConePos(2) = min(minConePos(2), min(datarun.cone_centers(datarun.cone_weights(:,ii)>0,2)));
            
            conesIncluded = max(conesIncluded, datarun.cone_weights(:,ii)>0); 
        end
    end
end

%plot(datarun.cone_centers(~~conesIncluded,1), datarun.cone_centers(~~conesIncluded,2), 'ko')

%% re-store selected RGCs in model structure

model.cone_centers = datarun.cone_centers;
model.rgc_COMs = datarun.rgc_COMs(includeNeuron,:);
model.cell_types = datarun.cell_types(includeNeuron);
model.cone_types = datarun.cone_types;
model.cone_weights = datarun.cone_weights(:,includeNeuron);
model.cones_included = conesIncluded;

% recenter x and y positions to be >= 0

offset = -minConePos + [100 100]; %add buffer region of 100 microns to avoid edge effects of spatial convolution with psf
model.cone_centers(:,1) = model.cone_centers(:,1) + ones(size(model.cone_centers(:,1)))*offset(1);
model.cone_centers(:,2) = model.cone_centers(:,2) + ones(size(model.cone_centers(:,1)))*offset(2);
model.rgc_COMs(:,1) = model.rgc_COMs(:,1) + ones(size(model.rgc_COMs(:,1)))*offset(1);
model.rgc_COMs(:,2) = model.rgc_COMs(:,2) + ones(size(model.rgc_COMs(:,1)))*offset(2);

%% plot spiders + zero s-cone weights
figure
hold on

for ii = 1:length(model.cell_types)
    center = model.rgc_COMs(ii,:);
    includeNeuron(ii) = true;
    plot(center(1), center(2), 'k.')
    coneWeightInds = find(model.cone_weights(:,ii));
    for jj = 1:length(coneWeightInds)        
        if model.cell_types(ii) == 3 %on midget
            plot([center(1) model.cone_centers(coneWeightInds(jj),1)],...
                [center(2) model.cone_centers(coneWeightInds(jj),2)], 'k-',...
                'LineWidth', abs(model.cone_weights(coneWeightInds(jj),ii))*2/max(model.cone_weights(:,ii)),...
                'Color', [0.5 0.5 0.5])
        elseif model.cell_types(ii) == 4 %off midget
            plot([center(1) model.cone_centers(coneWeightInds(jj),1)],...
                [center(2) model.cone_centers(coneWeightInds(jj),2)], 'k-',...
                'LineWidth', abs(model.cone_weights(coneWeightInds(jj),ii))*2/max(model.cone_weights(:,ii)),...
                'Color', [0 0 0])
        end
        
        if strcmpi(model.cone_types(coneWeightInds(jj)),'M')
            plot(model.cone_centers(coneWeightInds(jj),1),...
                model.cone_centers(coneWeightInds(jj),2), 'o', 'MarkerEdgeColor', [0.1 0.8 0.1])
        elseif strcmpi(model.cone_types(coneWeightInds(jj)),'L')
            plot(model.cone_centers(coneWeightInds(jj),1),...
                model.cone_centers(coneWeightInds(jj),2), 'o', 'MarkerEdgeColor', [1 0 0])
        else %if S cone, remove weight from model
            model.cone_weights(coneWeightInds(jj),ii) = 0;
            break
        end
        
        if ii == 22
            plot([center(1) model.cone_centers(coneWeightInds(jj),1)],...
                [center(2) model.cone_centers(coneWeightInds(jj),2)], 'k-',...
                'LineWidth', abs(model.cone_weights(coneWeightInds(jj),ii))*5/max(model.cone_weights(:,ii)),...
                'Color', [0.5 0.5 0.5])
        end
    end
end

%set(gca, 'xlim', [0 300], 'ylim', [0 250])
hold off
axis equal
ylabel('microns')
xlabel('microns')
set(gca, 'box', 'on')

%% generate stimulus

%in microns
spotCenter = [250, 200];
spotRadius = 30;
stimSize = [500 500];

%generates binary stimulus defining spatial locations of foreground (1) and background (0)
stimulus = generateSpotStimulus(spotCenter, spotRadius, stimSize);


%define spectra of foreground and background
%spectral phosphor power distributions: (:,1) = r, (:,2) = g, (:,3) = b
load B_monitor
phosphors = B_monitor;
spectrum = 380:5:780;

% make a achromatic stimulus: equal amount of r,g,b phosphors
gray_spec = phosphors(:,1) + phosphors(:,2) + phosphors(:,3);

% background/baseline stimulus = 1 x gray_spec, all stimuli relative to 1
bg_spec = gray_spec;

% foreground
fg_spec = gray_spec*3; %corresponds to 50% contrast


% stimulus time course (in ms): when 0, full-field backgro[und, when 1, fg
% applied
stimulus_tc = zeros(500,1); %must be multiple of 50
stimulus_tc(101:300) = 1;


figure
subplot(2,1,1)
imagesc(stimulus);
colormap(gray)
axis equal
set(gca, 'xlim', [0 500], 'ylim', [0 500])
title('stimulus foreground')
xlabel('microns')
ylabel('microns')

subplot(2,1,2)
plot(stimulus_tc,'k-')
xlabel('ms')
set(gca,'ylim', [-0.5 1.5], 'ytick', [0 1], 'yticklabel', {'bg alone','bg and fg'},'xlim', [0 500])
title('stimulus timecourse')

%% stimulus spatial components filtered by M or L psf

%separate foreground and background spatial regions
stimulus_fg_L_psf_filtered = psf_filter(stimulus, 'L');
stimulus_fg_M_psf_filtered = psf_filter(stimulus, 'M');
stimulus_bg_L_psf_filtered = psf_filter(~stimulus, 'L');
stimulus_bg_M_psf_filtered = psf_filter(~stimulus, 'M');

%full-field background
stimulus_ff_L_psf_filtered = psf_filter(ones(size(stimulus)), 'L');
stimulus_ff_M_psf_filtered = psf_filter(ones(size(stimulus)), 'M');


psfM = load('ap3mpsfSmall.txt');
psfL = load('ap3LpsfSmall.txt');

figure
subplot(3,2,1)
imagesc(psfL)
colormap(gray)
axis equal
set(gca, 'xlim', [0 100], 'ylim', [0 100])
title('L cone psf')
ylabel('microns')

subplot(3,2,2)
imagesc(psfM)
colormap(gray)
axis equal
set(gca, 'xlim', [0 100], 'ylim', [0 100])
title('M cone psf')

subplot(3,2,3)
imagesc(stimulus_fg_L_psf_filtered)
colormap(gray)
axis equal
set(gca, 'xlim', [0 500], 'ylim', [0 500])
title('foreground: L psf filtered')
ylabel('microns')

subplot(3,2,4)
imagesc(stimulus_fg_M_psf_filtered)
colormap(gray)
axis equal
set(gca, 'xlim', [0 500], 'ylim', [0 500])
title('foreground: M psf filtered')

subplot(3,2,5)
imagesc(stimulus_bg_L_psf_filtered)
colormap(gray)
axis equal
set(gca, 'xlim', [0 500], 'ylim', [0 500])
title('background: L psf filtered')
xlabel('microns')
ylabel('microns')

subplot(3,2,6)
imagesc(stimulus_bg_M_psf_filtered)
colormap(gray)
axis equal
set(gca, 'xlim', [0 500], 'ylim', [0 500])
title('background: M psf filtered')
xlabel('microns')

%% 
% determine difference in cone signal corresponding to 25% stimulus 
% contrast full-field stimulus (for baseline = 1x r+g+b spectra) for each
% cone type
%
% 25% stimulus contrast corresponds to 5/3x baseline stimulus 
%(Imax-Imin / Imax + Imin)

baseline_spec = phosphors(:,1) + phosphors(:,2) + phosphors(:,3);
max_spec = (5/3)*(phosphors(:,1) + phosphors(:,2) + phosphors(:,3));

ff_stimulus = ones(100,100);

L_filtered_stimulus = psf_filter(ff_stimulus, 'L');
M_filtered_stimulus = psf_filter(ff_stimulus, 'M');

L_cone_response_baseline = cone_filter([50 50], 'L', L_filtered_stimulus,...
    zeros(size(L_filtered_stimulus)), baseline_spec, zeros(size(baseline_spec)));
M_cone_response_baseline = cone_filter([50 50], 'M', M_filtered_stimulus,...
    zeros(size(M_filtered_stimulus)), baseline_spec, zeros(size(baseline_spec)));
L_cone_response_25 = cone_filter([50 50], 'L', L_filtered_stimulus,...
    zeros(size(L_filtered_stimulus)), max_spec, zeros(size(max_spec)));
M_cone_response_25 = cone_filter([50 50], 'M', M_filtered_stimulus,...
    zeros(size(M_filtered_stimulus)), max_spec, zeros(size(max_spec)));

L_cone_noise_std = L_cone_response_25 - L_cone_response_baseline;
M_cone_noise_std = M_cone_response_25 - M_cone_response_baseline;

%% generate noise vectors for participating cones
timeSegments = length(stimulus_tc)/50;
coneNoise = zeros(size(datarun.cone_centers,1), timeSegments);

includedConeInd = find(model.cones_included);


for ii = 1:length(includedConeInd)
    cone_type = model.cone_types(includedConeInd(ii));
    if strcmpi(cone_type,'L')
        coneNoise(includedConeInd(ii),:) = normrnd(0, L_cone_noise_std, 1, timeSegments);
    elseif strcmpi(cone_type,'M')
        coneNoise(includedConeInd(ii),:) = normrnd(0, M_cone_noise_std, 1, timeSegments);
    end
end

% plot a few cone noise timecourses
% figure
% hold on
% for ii = 1:10
%     cone_id = includedConeInd(ii);
%     temp = zeros(size(stimulus_tc));
%     for kk = 1:timeSegments
%         temp((kk-1)*50+1:kk*50) = coneNoise(cone_id,kk);
%     end
%     plot(temp + (ii-1)*0.05)
% end
% hold off
% title('example cone noise vectors')

%% filter stimulus for each cone and apply noise at 50 ms intervals
nCells = length(model.cell_types);

coneSumOutputClean = cell(nCells,1);
coneSumOutputNoisy = cell(nCells,1);
coneSumResponseProfile = cell(nCells,1);
coneSumResponseProfileNoWeights = cell(nCells,1);
coneSumResponseProfile_bl = cell(nCells,1);


%figure
%coneSumAxes = axes;

for ii = 1:nCells
%     ii = 20;
    coneSumOutputClean{ii} = zeros(size(stimulus_tc));
    coneSumOutputNoisy{ii} = zeros(size(stimulus_tc));
    coneSumResponseProfile{ii} = zeros(size(stimulus));
    coneSumResponseProfileNoWeights{ii} = zeros(size(stimulus));
    coneSumResponseProfile_bl{ii} = zeros(size(stimulus));
    coneWeightInds = find(model.cone_weights(:,ii));
    sumColors = hsv(length(coneWeightInds));
    for jj = 1:length(coneWeightInds)
        cone_pos = model.cone_centers(coneWeightInds(jj),:);
        cone_id = coneWeightInds(jj);
        cone_type = model.cone_types(cone_id);
        
        %during 'stimulus off' (full-field at background intensity)
        %apply point-spread function to stimulus
        
        %full-field background
        if strcmpi(cone_type, 'L')
            stimulus_ff = stimulus_ff_L_psf_filtered;
        elseif strcmpi(cone_type, 'M')
            stimulus_ff = stimulus_ff_M_psf_filtered;
        end
        
        
        %apply cone spatial and spectral filters
        [cone_response_baseline cone_response_bl_profile] = cone_filter(cone_pos, cone_type, zeros(size(stimulus_ff)),...
            stimulus_ff, fg_spec, bg_spec);
        
%         figure
%         imagesc(cone_response_bl_profile);
%         axis equal
        
        coneSumResponseProfile_bl{ii} = coneSumResponseProfile_bl{ii} + model.cone_weights(cone_id,ii)*cone_response_bl_profile;
        
        %during 'stimulus on'
        %apply point-spread function to stimulus
        if strcmpi(cone_type, 'L')
            stimulus_fg = stimulus_fg_L_psf_filtered;
            stimulus_bg = stimulus_bg_L_psf_filtered;
        elseif strcmpi(cone_type, 'M')
            stimulus_fg = stimulus_fg_M_psf_filtered;
            stimulus_bg = stimulus_bg_M_psf_filtered;
        end

        %apply cone spatial and spectral filters
        [cone_response_stimulus, cone_response_profile] = cone_filter(cone_pos, cone_type, stimulus_fg,...
            stimulus_bg, fg_spec, bg_spec);
        
%         figure
%         imagesc(cone_response_profile);
%         colormap(gray)
%         axis equal
%         title(['cone' num2str(cone_id)])
        
        coneSumResponseProfile{ii} = coneSumResponseProfile{ii} + model.cone_weights(cone_id,ii)*cone_response_profile;
        coneSumResponseProfileNoWeights{ii} = coneSumResponseProfileNoWeights{ii}+ cone_response_profile;
        
        % add cone noise and sum weighted cone outputs
        cone_response_tc = stimulus_tc*cone_response_stimulus +...
            (ones(size(stimulus_tc))-stimulus_tc)*cone_response_baseline;
        
        cone_noise_full = zeros(size(stimulus_tc));
        for kk = 1:timeSegments
            cone_noise_full((kk-1)*50+1:kk*50) = coneNoise(cone_id,kk);
        end
        
        cone_response_tc_noisy = cone_response_tc + cone_noise_full;
        
        coneSumOutputClean{ii} = coneSumOutputClean{ii} + model.cone_weights(cone_id,ii)*cone_response_tc;
        coneSumOutputNoisy{ii} = coneSumOutputNoisy{ii} + model.cone_weights(cone_id,ii)*cone_response_tc_noisy;

%         axes(coneSumAxes)
%         hold on
%         %plot(coneSumOutputNoisy{ii}, 'color', sumColors(jj,:))
% %         plot(model.cone_weights(cone_id,ii)*cone_response_tc_noisy, 'color', sumColors(jj,:))
%         plot(model.cone_weights(cone_id,ii)*cone_response_tc, 'color', sumColors(jj,:))
%         %plot(coneSumOutputClean{ii}, 'color', sumColors(jj,:), 'LineWidth', 2)
%         hold off
    end
%     axes(coneSumAxes)
%     hold on
% %     plot(coneSumOutputNoisy{ii}, 'k')
%     plot(coneSumOutputClean{ii}, 'k')
%     hold off
%     xlabel('time (ms)')
%     ylabel('cone response (a.u.)')
%     
    
%     figure
%     hold on
%     imagesc(coneSumResponseProfile{ii},[0 15])
%     %imagesc(coneSumResponseProfile{ii} - coneSumResponseProfile_bl{ii})
%     %imagesc(coneSumResponseProfileNoWeights{ii})
%     colormap(gray)
%     %superimpose original stimulus
%     th = 0:0.01:2*pi;
%     xCircle = cos(th)*spotRadius + spotCenter(1);
%     yCircle = sin(th)*spotRadius + spotCenter(2);
%     plot(yCircle, xCircle, 'w-')
%     hold off
%     axis equal
%     title(['cell' num2str(ii)])    

end

% maxResponse = 0;
% allCellsConeSumResponseProfile = zeros(size(coneSumResponseProfile{1}));
% for ii = 1:nCells
%     allCellsConeSumResponseProfile = allCellsConeSumResponseProfile + coneSumResponseProfile{ii};
%     %figure
%     %imagesc(coneSumResponseProfile{ii}, [0 1])
%     %title(num2str(ii))
%     maxResponse = max(maxResponse, max(max(coneSumResponseProfile{ii})));
% end

% th = 0:0.01:2*pi;
% xCircle = cos(th)*spotRadius + spotCenter(1);
% yCircle = sin(th)*spotRadius + spotCenter(2);
% 

% figure
% hold on
% imagesc(allCellsConeSumResponseProfile, [0 50])
% plot(yCircle, xCircle, 'w-')
% hold off
% axis equal


%%
return

%% RF time course
%from glm tutorial
% Create a default (temporal) stimulus filter

nkt = 20;
tk = [0:nkt-1]';
b1 = nkt/32; b2 = nkt/16;
k1 = 1/(gamma(6)*b1)*(tk/b1).^5 .* exp(-tk/b1);  % Gamma pdfn
k2 = 1/(gamma(6)*b2)*(tk/b2).^5 .* exp(-tk/b2);  % Gamma pdf
tc_filter_forward = k1-k2./1.5;
%resample at each ms (peak ~30 ms)
tc_filter_forward = interp(tc_filter_forward, 10);

tc_filter = flipud(tc_filter_forward);

% temporally filter sum of cone signals and apply nonlinearity
t_filtered_response = cell(nCells,1);
min_response = 10000;
max_response = 0;
for ii = 1:nCells
    %t_filtered_response{ii} = conv(coneSumOutputNoisy{ii}, tc_filter_forward, 'same');
    t_filtered_response{ii} = sameconv(coneSumOutputNoisy{ii}, tc_filter);
    min_response = min(min_response, min(t_filtered_response{ii}));
    max_response = max(max_response, max(t_filtered_response{ii}));
end

% scale responses so that they fall in the range from 2 to 4 and resample
% at 0.1 ms
figure
hold on
StimHi = cell(nCells,1);
for ii = 1:nCells
    t_filtered_response{ii} = (t_filtered_response{ii}-min_response)*(8/max_response);
    StimHi{ii} = interp(t_filtered_response{ii}(1:500), 10);
    StimHi{ii}(1:800) = 0;
    plot(StimHi{ii})
end
hold off

%% from glm tutorial

ggsim.nlfun = @exp;
DTsim = 0.01;
RefreshRate = 100;

% post-spike current (from glm tutorial)
dt = 0.01; %simulation time steps

%version 1
% --- Make basis for post-spike (h) current ------
% ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
% ihbasprs.hpeaks = [.1 2];  % Peak location for first and last vectors
% ihbasprs.b = .5;  % How nonlinear to make spacings
% ihbasprs.absref = .1; % absolute refractory period 
% [iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs, dt);
% ih = ihbasis*[-10 -5 0 2 -2]';  % h current



% Make basis for (h) 
ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [0 3];  % Peak location for first and last vectors
ihbasprs.b = .5;  % How nonlinear to make spacings
ihbasprs.absref = []; % absolute refractory period (optional param)
[iht,ihbasOrthog,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim);

% store these in the simulation structure
ggsim.ihbasprs = ihbasprs;  
ggsim.iht = iht;
% Make post-spike filter h using this basis
ihweights = [-10 1 1 -1 -.5]';  % weights
ih = ihbasis*ihweights;  % post-spike filter


% Step 5: finally, we need a loop to generate spikes and add in the
% post-spike filter contribution every time there's a spike
nTimeBins = length(StimHi{1}); % total length of high-frequency time bins
hlen = length(ih);% length of post-spike filter

stim_with_filter = cell(nCells,1);
sptimes = cell(nCells,1);
nlValues = cell(nCells,1);
%for ii = 1:nCells
ii = 20;
    stim_with_filter{ii} = StimHi{ii};  % This stores the total filter output (stim + spike filter)
    sptimes{ii} = [];
    nlValues{ii} = zeros(size(StimHi{ii}));

    % Loop: apply nonlinearity to net filter output and draw spikes
    for j = 1:nTimeBins
        % Draw a bernoulli random variable to determine whether to spike or not
        nlValues{ii}(j) = ggsim.nlfun(stim_with_filter{ii}(j));
        if rand < ggsim.nlfun(stim_with_filter{ii}(j))*DTsim/RefreshRate %probability of spiking = nlfun value * dt (in s)
            % Spike!
            sptimes{ii} = [sptimes{ii}; DTsim*j];
            % add h to subsequent time bins
        stim_with_filter{ii}(j+1:min(nTimeBins,j+hlen)) = stim_with_filter{ii}(j+1:min(nTimeBins,j+hlen)) + ...
                ih(1:length(j+1:min(nTimeBins,j+hlen)));
        end
    end
%end

% That's all!  Plot it:
tthi = (1:nTimeBins)/10;

subplot(311)
plot(tthi)
plot(1:tthi(end), coneSumOutputNoisy{ii}(1:500))
title('summed cone signal')
ylabel('a.u.')

subplot(312);
plot(tthi, StimHi{ii})
% plot(tthi, StimHi{ii}, tthi, stim_with_filter{ii});
title('temporal filter output');
axis tight;
subplot(313);
plot(tthi, stim_with_filter{ii}-StimHi{ii}, sptimes{ii}*10, 1, 'r.'); %just plots results of post-spike filter and spikes
title('time filter output and spikes');
axis tight;
xlabel('time (ms)')
