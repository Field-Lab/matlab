% Attempt to simulate the C/R response of midget versus parasol RGCs based
% on C/R of single cones and STA/RF weights.

% Generally seems to give higher midget RGC C/R than expected based on
% Kaplan and Shapley parvocellular results.


%% Basics
piece = '2012-09-13-2';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
keep_vars = {'piece'; 'loadopts'; 'keep_vars'; 'd'};

%% Streamed runs and cell picks
streamedopts = loadopts;
streamedopts.load_sta = true;
streamedopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
streamedopts.set_polarities = {'guess', true};

d.d05s = load_data([piece '/streamed/data005/data005'], streamedopts);
d.d05s.offM = [1727 1816 1966 2026 2116 2386 2942 6991];
d.d05s.offP = [4441 5371 5551 6841];

d.d05s = load_cones(d.d05s);

%% Load 1cone stim run
d.d07 = load_data([piece '/data007'], loadopts);

d.d07 = read_stim_lisp_output(d.d07);
d.d07.stimulus = parse_stim_rgbs(d.d07.stimulus);
d.d07.stimulus = parse_cr_rgbs(d.d07.stimulus);

cell_list_map = map_ei(d.d05s, d.d07, 'master_cell_type', 'all');
d.d07.offM_fd05s = cell2mat(cell_list_map(get_cell_indices(d.d05s, d.d05s.offM)));
d.d07.offP_fd05s = cell2mat(cell_list_map(get_cell_indices(d.d05s, d.d05s.offP)));
clear cell_list_map;

%% Setup
conerun = d.d05s;
conerun.rgcs = [conerun.offM conerun.offP];
rasterrun = d.d07;
rasterrun.rgcs = [rasterrun.offM_fd05s rasterrun.offP_fd05s];
nM = length(conerun.offM);
nP = length(conerun.offP);
nrgcs = nM + nP;

%% Boxcar template
boxstart = 0.05;
boxend = 0.25;
template.x = -0.1:0.02:0.5;
template.y = zeros(size(template.x));
template.y(template.x > boxstart & template.x < boxend) = 1;
template.y = template.y ./ sum(template.y);

%% Calc the whole block of C/R values, rearrange into convenient dimensionality
stimulus = rasterrun.stimulus;
[crs rasterhistsy] = single_map_contrast_response(rasterrun, rasterrun.triggers(1:2:end), rasterrun.rgcs, template, cell2mat(stimulus.cr.single_cones'));
crs = permute(reshape(crs, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
rasterhistsy = permute(reshape(rasterhistsy, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
crsx = cell2mat(rasterrun.stimulus.cr.single_cone_intensities);
ncontrasts = size(crs,2); % Assume same for all

%% Pick cone to use just based on max response
for k = 1:nrgcs
    [~,maxcones(k)] = max(sum(crs(:,:,k),2));
    [~,mincones(k)] = min(sum(crs(:,:,k),2));
end

% Manual selection
mancones = maxcones;
mancones(10) = 1;

%% Fit max/min/man cone with NORMCDF
x = crsx(1,:)./3;

fitfunc = @(p,x) (p(4) + p(3).*normcdf(x, p(1), p(2)));
p0 = [0.5 0.25 40 5];
optimopts = optimset('Algorithm', 'levenberg-marquardt', 'MaxFunEvals', 5000, 'MaxIter', 5000);
for k = 1:nrgcs
    maxcone = maxcones(k);
    y = crs(maxcone,:,k);    
    psmax{k} = lsqcurvefit(@(p,x)(p(4) + p(3).*normcdf(x,p(1),p(2))), p0, -x, y, [], [], optimopts);
    
    y = crs(mincones(k),:,k);    
    psmin{k} = lsqcurvefit(@(p,x)(p(4) + p(3).*normcdf(x,p(1),p(2))), p0, -x, y, [], [], optimopts);
    
    y = crs(mancones(k),:,k);    
    psman{k} = lsqcurvefit(@(p,x)(p(4) + p(3).*normcdf(x,p(1),p(2))), p0, -x, y, [], [], optimopts);
end

% Plot?
if true
    figure();
    hold on
    
    fitstep = 0.01;
    fitx = x(1):fitstep:x(end);
    for k = 1:nM
        y = crs(mancones(k),:,k);
        p = psman{k};
        
        plot(x, y, 'ob');
        fit = p(4) + p(3).*normcdf(-fitx, p(1), p(2));
        plot(fitx, fit, 'b');
    end
    for k = nM + (1:nP)
        y = crs(mancones(k),:,k);
        p = psman{k};
        
        plot(x, y, 'sr');
        fit = p(4) + p(3).*normcdf(-fitx, p(1), p(2));
        plot(fitx, fit, '--r');
    end
end

%% Get cone weights for center
coneweightopts = struct('thresh', 0.05, 'radius', [0 10], 'polarity', 1, 'contiguity', false,'scale', 3.0);
[mosaic_weights, selection, extras] = select_cone_weights(conerun, conerun.rgcs, coneweightopts);

%% Figure out which cone in each stimulus index is actually intended for
% each RGC.  Multiple ways to do this but mosaic_weights is good enough.
cones_stimulated = recover_cones_stimulated(conerun, rasterrun.stimulus.mapims{1});
ncones = length(cones_stimulated);
for i = 1:ncones
    cones = cones_stimulated{i};
    coneweights = mosaic_weights(cones,:);
    [~,maxinds] = max(coneweights);
    rasterrun.stimconenums(:,i) = cones_stimulated{i}(maxinds);
end

%% Check fits scaled by cone weight; not technically correct scaling as we
% actually should use the integral of RF within the region, but this should
% be reasonably close
% colors = jet(ncones);
% fitstep = 0.01;
% fitx = x(1):fitstep:x(end);
% for k = 1:nrgcs
%     figure;
%     hold on
% 
%     stimcones = rasterrun.stimconenums(k,:);
%     weights = mosaic_weights(stimcones,k);
% 
%     maxcone = maxcones(k);
%     weights = weights ./ weights(maxcone);
%     for i = 1:ncones
%         y = crs(i,:,k);
%         
%         h = plot(x, y, 'o');
%         color = colors(i,:);
%         set(h, 'Color', color);
%         
%         p = ps{k};
%         fit = p(4) + p(3).*normcdf(-fitx.*weights(i), p(1), p(2));
%         h = plot(fitx, fit);
%         set(h, 'Color', color);
%     end
% end
% Not so close...  But let's move on.

%% Simulate C/R for whole RGC, assuming C/R of individual cones match C/R fit and cones sum linearly
normcones = mancones;
ps = psman;

x0 = find(crsx(1,:) == 0);

simstep = 0.01;
simx = x(1):simstep:x(end);
numx = length(simx);
for k = 1:nrgcs
    coneweights = mosaic_weights(selection(:,k),k);
    normconeweight = mosaic_weights(rasterrun.stimconenums(k,normcones(k)), k); % Should match max(coneweights)...
    normconeweights = coneweights ./ normconeweight;
    ncones = length(normconeweights);
    
    simxes = normconeweights * simx;
    p = ps{k};
    simys = p(4) + p(3).*normcdf(-simxes, p(1), p(2));
    baseline = p(4) + p(3).*normcdf(0, p(1), p(2));
    simy{k} = sum(simys) - (ncones-1).*baseline;
end

offMCRs = cell2mat(simy(1:8)');
offPCRs = cell2mat(simy(9:12)');
figure; hold on
plot(simx, offMCRs(:,:), 'b');
plot(simx, offPCRs(:,:), '--r');

%% Simulate C/R for whole RGC, assuming C/R of individual cones match C/R fit and cones sum linearly
normcones = mancones;
ps = psman;

x0 = find(crsx(1,:) == 0);

simstep = 0.01;
simx = x(1):simstep:x(end);
numx = length(simx);
for k = 1:nrgcs
    coneweights = mosaic_weights(selection(:,k),k);
    normconeweight = mosaic_weights(rasterrun.stimconenums(k,normcones(k)), k);
    normconeweights = coneweights ./ normconeweight;
    
    scaleweights{k} = sum(normconeweights);
    p = ps{k};
    simy{k} = p(4) + p(3).*normcdf(-simx .* scaleweights{k}, p(1), p(2));
end

offMCRs = cell2mat(simy(1:8)');
offPCRs = cell2mat(simy(9:12)');
figure; hold on
plot(simx, offMCRs(:,:), 'b');
plot(simx, offPCRs(:,:), '--r');

for k = 1:nrgcs
    if k <= nM, style = 'ob'; else style = 'sr'; end
    plot(crsx(1,:)./scaleweights{k}, crs(normcones(k),:,k), style);
end

%% Just data
for k = 1:nrgcs
    for i = 1:4
        coneweights = mosaic_weights(selection(:,k),k);
        normconeweight = mosaic_weights(rasterrun.stimconenums(k,i), k);
        normconeweights = coneweights ./ normconeweight;
        scaleweights(k,i) = sum(normconeweights);
    end
end

figure; hold on
for k = 1:nrgcs
    if k <= nM, style = '-b'; else style = '--r'; end
    for i = 4 % 1:4
        y = crs(i,:,k);
%        y = y - y(5);
        plot(crsx(1,:)./scaleweights(k,i), y, style);
    end
end