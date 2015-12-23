% Attempt to simulate the C/R response of midget versus parasol RGCs based
% on C/R of single cones and STA/RF weights.

% Simpler version using only real data not fits, as suggested by EJ.

%% Basics
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
streamedopts = loadopts;
streamedopts.load_sta = true;
streamedopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
streamedopts.set_polarities = {'guess', true};
pieces = containers.Map('KeyType', 'char', 'ValueType', 'any');
keep_vars = {'pieces'; 'loadopts'; 'streamedopts'; 'keep_vars'};

%% 2011-06-24-6
piece = '2011-06-24-6';
d.d05s = load_data([piece '/streamed/data005/data005'], streamedopts);
d.d08s = load_data([piece '/streamed/data008/data008'], streamedopts);
d.d09  = load_data([piece '/data009'], loadopts);
d.d11  = load_data([piece '/data011'], loadopts);
d.d13  = load_data([piece '/data013'], streamedopts);
d.d09.mapd05s = map_ei(d.d05s, d.d09);
d.d13.mapd05s = map_ei(d.d05s, d.d13);
d.d11.mapd08s = map_ei(d.d08s, d.d11);
d.d13.mapd08s = map_ei(d.d08s, d.d13);

% Convert old style additivity maps into format expected by cr_compound_plot
% Map 0 and 1 are the two cones, map 2 is both
d.d09 = read_stim_lisp_output(d.d09, '2011-06-24-6_data005_1and2');
    stimstruct = d.d09.stimulus;
    stimstruct.numcones = 2;
    stimstruct.mapims = {combine_maps(stimstruct.mapims(1:2))};
    stimstruct.mapnyc     = {cellfun(@(C)(C{1}), stimstruct.mapnyc(1:2),     'UniformOutput', false)};
    stimstruct.mapnycpoly = {cellfun(@(C)(C{1}), stimstruct.mapnycpoly(1:2), 'UniformOutput', false)};
    rgbs = stimstruct.rgbs;
    stimstruct.rgbs = cell(1, length(stimstruct.pulses));
    stimstruct.rgbs(stimstruct.maps == 0) = cellfun(@(A)([A; 0 0 0]), mat2cell(rgbs(stimstruct.maps == 0,:), ones(1,sum(stimstruct.maps == 0)), 3), 'UniformOutput', false);
    stimstruct.rgbs(stimstruct.maps == 1) = cellfun(@(A)([0 0 0; A]), mat2cell(rgbs(stimstruct.maps == 1,:), ones(1,sum(stimstruct.maps == 1)), 3), 'UniformOutput', false);
    stimstruct.rgbs(stimstruct.maps == 2) = cellfun(@(A)([A;     A]), mat2cell(rgbs(stimstruct.maps == 2,:), ones(1,sum(stimstruct.maps == 2)), 3), 'UniformOutput', false);
    d.d09.stimulus = stimstruct;
        d.d11 = read_stim_lisp_output(d.d11, '2011-06-24-6_f08_1and2');
        stimstruct = d.d11.stimulus;
        stimstruct.numcones = 2;
        stimstruct.mapims = {combine_maps(stimstruct.mapims(1:2))};
        stimstruct.mapnyc     = {cellfun(@(C)(C{1}), stimstruct.mapnyc(1:2),     'UniformOutput', false)};
        stimstruct.mapnycpoly = {cellfun(@(C)(C{1}), stimstruct.mapnycpoly(1:2), 'UniformOutput', false)};
        rgbs = stimstruct.rgbs;
        stimstruct.rgbs = cell(1, length(stimstruct.pulses));
        stimstruct.rgbs(stimstruct.maps == 0) = cellfun(@(A)([A; 0 0 0]), mat2cell(rgbs(stimstruct.maps == 0,:), ones(1,sum(stimstruct.maps == 0)), 3), 'UniformOutput', false);
        stimstruct.rgbs(stimstruct.maps == 1) = cellfun(@(A)([0 0 0; A]), mat2cell(rgbs(stimstruct.maps == 1,:), ones(1,sum(stimstruct.maps == 1)), 3), 'UniformOutput', false);
        stimstruct.rgbs(stimstruct.maps == 2) = cellfun(@(A)([A;     A]), mat2cell(rgbs(stimstruct.maps == 2,:), ones(1,sum(stimstruct.maps == 2)), 3), 'UniformOutput', false);
        d.d11.stimulus = stimstruct;
d.d09.stimulus = parse_stim_rgbs(d.d09.stimulus);
d.d11.stimulus = parse_stim_rgbs(d.d11.stimulus);
d.d09.stimulus = parse_cr_rgbs(d.d09.stimulus);
d.d11.stimulus = parse_cr_rgbs(d.d11.stimulus);
pieces(piece) = d;
leave(keep_vars{:});

piece = '2011-06-24-6';
d = pieces('2011-06-24-6');
midgetconerun = d.d05s;
midgetconerun.rgcs = [50 78 64 137 158 120 15 29 42];
parasolconerun = d.d08s;
parasolconerun.rgcs = [34 75 88 96 125 198];
midgetrun = d.d09;
midgetrun.rgcs = d.d09.mapd05s(get_cell_indices(midgetconerun, midgetconerun.rgcs));
parasolrun = d.d11;
parasolrun.rgcs = d.d11.mapd08s(get_cell_indices(parasolconerun, parasolconerun.rgcs));

% Boxcar template
boxstart = 0.05;
boxend = 0.25;
template.x = -0.1:0.02:0.5;
template.y = zeros(size(template.x));
template.y(template.x > boxstart & template.x < boxend) = 1;
template.y = template.y ./ sum(template.y);

% Calc the whole block of C/R values, rearrange into convenient dimensionality
rasterrun = midgetrun;
stimulus = rasterrun.stimulus;
[crs rasterhistsy] = single_map_contrast_response(rasterrun, rasterrun.triggers, rasterrun.rgcs, template, cell2mat(stimulus.cr.single_cones'));
crs = permute(reshape(crs, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
rasterhistsy = permute(reshape(rasterhistsy, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
crsx = cell2mat(rasterrun.stimulus.cr.single_cone_intensities);
ncontrasts = size(crs,2); % Assume same for all
midgetcrs = crs;
midgethistsy = rasterhistsy;

% Repeat for parasols
rasterrun = parasolrun;
stimulus = rasterrun.stimulus;
[crs rasterhistsy] = single_map_contrast_response(rasterrun, rasterrun.triggers, rasterrun.rgcs, template, cell2mat(stimulus.cr.single_cones'));
crs = permute(reshape(crs, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
rasterhistsy = permute(reshape(rasterhistsy, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
parasolcrs = crs;
parasolhistsy = rasterhistsy;

% Get cone weights for centers
coneweightopts = struct('thresh', 0.05, 'radius', [0 10], 'polarity', 1, 'contiguity', false,'scale', 3.0);
midgetconerun  = load_cones(midgetconerun, 1);
parasolconerun = load_cones(parasolconerun);
[midget_mosaic_weights,  midget_selection]  = select_cone_weights(midgetconerun,  midgetconerun.rgcs,  coneweightopts);
[parasol_mosaic_weights, parasol_selection] = select_cone_weights(parasolconerun, parasolconerun.rgcs, coneweightopts);


% Figure out which cone in each stimulus index is actually intended for
% each RGC.  Multiple ways to do this but mosaic_weights is good enough.
midget_cones_stimulated  = recover_cones_stimulated(midgetconerun,  midgetrun.stimulus.mapims{1});
parasol_cones_stimulated = recover_cones_stimulated(parasolconerun, parasolrun.stimulus.mapims{1});
for i = 1:length(midget_cones_stimulated);
    cones = midget_cones_stimulated{i};
    coneweights = midget_mosaic_weights(cones,:);
    [~,maxinds] = max(coneweights);
    midgetrun.stimconenums(:,i) = midget_cones_stimulated{i}(maxinds);
end
for i = 1:length(parasol_cones_stimulated);
    cones = parasol_cones_stimulated{i};
    coneweights = parasol_mosaic_weights(cones,:);
    [~,maxinds] = max(coneweights);
    parasolrun.stimconenums(:,i) = parasol_cones_stimulated{i}(maxinds);
end

% Calculate scaling based on mosaic_weight
for k = 1:length(midgetconerun.rgcs);
    coneweights = midget_mosaic_weights(midget_selection(:,k),k);
    normconeweight = midget_mosaic_weights(midgetrun.stimconenums(k,1), k);
    normconeweights = coneweights ./ normconeweight;
    midgetscaleweight(k) = sum(normconeweights);
end
for k = 1:length(parasolconerun.rgcs);
    coneweights = parasol_mosaic_weights(parasol_selection(:,k),k);
    normconeweight = parasol_mosaic_weights(parasolrun.stimconenums(k,1), k);
    normconeweights = coneweights ./ normconeweight;
    parasolscaleweight(k) = sum(normconeweights);
end

% Plot
figure
subplot(1, 2, 1); hold on
for k = 1:length(midgetconerun.rgcs);
    y = midgetcrs(1,:,k);
%    y = y - y(5);
    plot(crsx(1,:), y, '-b');
end
for k = 1:length(parasolconerun.rgcs);
    y = parasolcrs(1,:,k);
%    y = y - y(5);
    plot(crsx(1,:), y, '--r');
end
subplot(1, 2, 2); hold on
for k = 1:length(midgetconerun.rgcs);
    y = midgetcrs(1,:,k);
%    y = y - y(5);
    plot(crsx(1,:)./midgetscaleweight(k), y, '-b');
end
for k = 1:length(parasolconerun.rgcs);
    y = parasolcrs(1,:,k);
%    y = y - y(5);
    plot(crsx(1,:)./parasolscaleweight(k), y, '--r');
end

%% 2012-09-13-2
piece = '2012-09-13-2';
conerun = load_data([piece '/streamed/data005/data005'], streamedopts);
conerun.offM = [1727 1816 1966 2026 2116 2386 2942 6991];
conerun.offP = [4441 5371 5551 6841];
conerun.rgcs = [conerun.offM conerun.offP];
conerun.conepicks = [4 4 4 4 4 4 4 4 4 4 4 4];
rasterrun = load_data([piece '/data007'], loadopts);
rasterrun.map = map_ei(conerun, rasterrun);
rasterrun.offM = rasterrun.map(get_cell_indices(conerun, conerun.offM));
rasterrun.offP = rasterrun.map(get_cell_indices(conerun, conerun.offP));
rasterrun.rgcs = [rasterrun.offM rasterrun.offP];
rasterrun = read_stim_lisp_output(rasterrun);
rasterrun.stimulus = parse_stim_rgbs(rasterrun.stimulus);
rasterrun.stimulus = parse_cr_rgbs(rasterrun.stimulus);
pieces('piece') = struct('conerun', conerun, 'rasterrun', rasterrun);

% Boxcar template
boxstart = 0.05;
boxend = 0.25;
template.x = -0.1:0.02:0.5;
template.y = zeros(size(template.x));
template.y(template.x > boxstart & template.x < boxend) = 1;
template.y = template.y ./ sum(template.y);

% Calc the whole block of C/R values, rearrange into convenient dimensionality
stimulus = rasterrun.stimulus;
[crs rasterhistsy] = single_map_contrast_response(rasterrun, rasterrun.triggers(1:2:end), rasterrun.rgcs, template, cell2mat(stimulus.cr.single_cones'));
crs = permute(reshape(crs, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
rasterhistsy = permute(reshape(rasterhistsy, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
crsx = cell2mat(rasterrun.stimulus.cr.single_cone_intensities);
ncontrasts = size(crs,2); % Assume same for all

% Get cone weights for center
conerun = load_cones(conerun);
coneweightopts = struct('thresh', 0.05, 'radius', [0 10], 'polarity', 1, 'contiguity', false,'scale', 3.0);
[mosaic_weights, selection, extras] = select_cone_weights(conerun, conerun.rgcs, coneweightopts);

% Figure out which cone in each stimulus index is actually intended for
% each RGC.  Multiple ways to do this but mosaic_weights is good enough.
cones_stimulated = recover_cones_stimulated(conerun, rasterrun.stimulus.mapims{1});
ncones = length(cones_stimulated);
for i = 1:ncones
    cones = cones_stimulated{i};
    coneweights = mosaic_weights(cones,:);
    [~,maxinds] = max(coneweights);
    rasterrun.stimconenums(:,i) = cones_stimulated{i}(maxinds);
end

% Calculate scaling based on mosaic_weight
nrgcs = length(conerun.rgcs);
for k = 1:nrgcs
    coneweights = mosaic_weights(selection(:,k),k);
    conepick = conerun.conepicks(k);
    normconeweight = mosaic_weights(rasterrun.stimconenums(k,conepick), k);
    normconeweights = coneweights ./ normconeweight;
    scaleweight(k) = sum(normconeweights);
end

% Plot
nM = length(conerun.offM);
nP = length(conerun.offP);
figure
subplot(1, 2, 1); hold on
for k = 1:nrgcs
    conepick = conerun.conepicks(k);
    if k <= nM, style = '-b'; else style = '--r'; end
    y = crs(conepick,:,k);
%    y = y - y(5);
    plot(crsx(1,:), y, style);
end
subplot(1, 2, 2); hold on
for k = 1:nrgcs
    conepick = conerun.conepicks(k);
    if k <= nM, style = '-b'; else style = '--r'; end
    y = crs(conepick,:,k);
%    y = y - y(5);
    plot(crsx(1,:)./scaleweight(k), y, style);
end

%% 2012-09-24-5
% Not that clean, only 2 negative contrasts

% piece = '2012-09-24-5';
% d.d03s = load_data([piece '/data003'], streamedopts);
% d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
% d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
% d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], streamedopts);
% d.d06cm = read_stim_lisp_output(d.d06cm);
% d.d03_06_07.mapd03s = map_ei(d.d03s, d.d03_06_07);
% pieces(piece) = d;
% leave(keep_vars{:});
% 
% piece = '2012-09-24-5';
% d = pieces('2012-09-24-5');
% conerun = d.d03s;
% conerun.offM = [3677 4591 3631 3646 5206 6931 7231 856];
% conerun.offP = [4172 3647 887 722 2971 3616 2521];
% conerun.rgcs = [conerun.offM conerun.offP];
% conerun.conepicks = [3 3 4 3 4 2 1 4        2 4 3 1 1 1 4];
% rasterrun = d.d06cm;
% rasterrun.rgcs = d.d03_06_07.mapd03s(get_cell_indices(conerun, conerun.rgcs));
% rasterrun = read_stim_lisp_output(rasterrun);
% rasterrun.stimulus = parse_stim_rgbs(rasterrun.stimulus);
% rasterrun.stimulus = parse_cr_rgbs(rasterrun.stimulus);
% 
% % Boxcar template
% boxstart = 0.05;
% boxend = 0.25;
% template.x = -0.1:0.02:0.5;
% template.y = zeros(size(template.x));
% template.y(template.x > boxstart & template.x < boxend) = 1;
% template.y = template.y ./ sum(template.y);
% 
% % Calc the whole block of C/R values, rearrange into convenient dimensionality
% stimulus = rasterrun.stimulus;
% [crs rasterhistsy] = single_map_contrast_response(rasterrun, rasterrun.triggers(1:2:end), rasterrun.rgcs, template, cell2mat(stimulus.cr.single_cones'));
% crs = permute(reshape(crs, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
% rasterhistsy = permute(reshape(rasterhistsy, length(rasterrun.rgcs), [], stimulus.numcones), [3 2 1]);
% crsx = cell2mat(rasterrun.stimulus.cr.single_cone_intensities);
% ncontrasts = size(crs,2); % Assume same for all
% 
% % Get cone weights for center
% conerun = load_cones(conerun);
% coneweightopts = struct('thresh', 0.05, 'radius', [0 10], 'polarity', 1, 'contiguity', false,'scale', 3.0);
% [mosaic_weights, selection, extras] = select_cone_weights(conerun, conerun.rgcs, coneweightopts);
% 
% % Figure out which cone in each stimulus index is actually intended for
% % each RGC.  Multiple ways to do this but mosaic_weights is good enough.
% cones_stimulated = recover_cones_stimulated(conerun, rasterrun.stimulus.mapims{1});
% ncones = length(cones_stimulated);
% for i = 1:ncones
%     cones = cones_stimulated{i};
%     coneweights = mosaic_weights(cones,:);
%     [~,maxinds] = max(coneweights);
%     rasterrun.stimconenums(:,i) = cones_stimulated{i}(maxinds);
% end
% 
% % Calculate scaling based on mosaic_weight
% nrgcs = length(conerun.rgcs);
% for k = 1:nrgcs
%     coneweights = mosaic_weights(selection(:,k),k);
%     conepick = conerun.conepicks(k);
%     normconeweight = mosaic_weights(rasterrun.stimconenums(k,conepick), k);
%     normconeweights = coneweights ./ normconeweight;
%     scaleweight(k) = sum(normconeweights);
% end
% 
% % Plot
% nM = length(conerun.offM);
% nP = length(conerun.offP);
% figure
% subplot(1, 2, 1); hold on
% for k = 1:nrgcs
%     conepick = conerun.conepicks(k);
%     if k <= nM, style = '-b'; else style = '--r'; end
%     y = crs(conepick,:,k);
% %    y = y - y(5);
%     plot(crsx(1,:), y, style);
% end
% subplot(1, 2, 2); hold on
% for k = 1:nrgcs
%     conepick = conerun.conepicks(k);
%     if k <= nM, style = '-b'; else style = '--r'; end
%     y = crs(conepick,:,k);
% %    y = y - y(5);
%     plot(crsx(1,:)./scaleweight(k), y, style);
% end
