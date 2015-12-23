% Basics
piece = '2012-09-24-5';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
keep_vars = {'piece'; 'loadopts'; 'staopts'; 'keep_vars'; 'd'};


%% Nice EIs
d.d01 = load_data([piece '/data001'], staopts);
d.d01.offM = [2431 2492 2641 1277 1338 1531 2252 2356 2386 2401];
d.d01.offP = [1324 2042 2521];


%% Plot EIs
% for i = 1:length(d.d01.offM)
%     figure();
%     subplot(2,1,1);
%     plot_ei_(get_ei(d.d01, d.d01.offM(i)), d.d01.ei.position);
%     
%     subplot(2,1,2);
%     plot_ei_(get_csd(d.d01, d.d01.offM(i)), d.d01.ei.position);
% end

for i = 1:length(d.d01.offP)
    figure();
    subplot(2,1,1);
    plot_ei_(get_ei(d.d01, d.d01.offP(i)), d.d01.ei.position);
    
    subplot(2,1,2);
    plot_ei_(get_csd(d.d01, d.d01.offP(i)), d.d01.ei.position);
end


%% EI figure
i = 2;
points = 1001;
offset = 8e-3;
elecs = [141 118 107 462 406 278];
startspike = 65;

cid = d.d01.offP(i);
cnum = get_cell_indices(d.d01, cid);
spikes = d.d01.spikes{cnum};

% Helper plots
% figure();
% semilogy(diff(spikes), '.-');
% hold on;
% semilogy(1:length(spikes), points / d.d01.sampling_rate, '-k');


% Get raw data
rawpath = fullfile(server_data_path, 'rawdata', piece, 'data001.bin');
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(rawpath);
start = (spikes(startspike) - offset) * d.d01.sampling_rate;
rawdata = rawFile.getData(start, points);
rawFile.close();
rawt = ((1:points) - 1)/d.d01.sampling_rate * 1e3;

% Manufactured spike trace
localspiketimes = spikes - spikes(startspike) + offset;
localspikes = find(localspiketimes > rawt(1)/1e3 & localspiketimes < rawt(end)/1e3);
localspikeidx = round(localspiketimes(localspikes)*20000);
spiket = ((1:0.5:points) - 1)/d.d01.sampling_rate * 1e3;
spiketrace = zeros(size(spiket));
spiketrace(localspikeidx*2) = 1;
spiketrace(localspikeidx*2-1) = 1;
spiketrace(localspikeidx*2+1) = 1;



% Plot...
clf();
set(gcf(), 'Color', 'w');

sanesubplot(5+length(elecs)*5, 8, {5 1:3});
plot(spiket, spiketrace, 'k');
axis off

ei = get_ei(d.d01, cid);
eiidx = 1:50;
eit = (eiidx - 1 - d.d01.ei.nlPoints)/d.d01.sampling_rate * 1e3;

% Plot raw traces
ylim = [-300 250];
for i = 1:length(elecs)
    rawax(i) = sanesubplot(1+length(elecs), 8, {1+i 1:3});
    
    % Highlights
    for j = 1:length(localspikes)
        st = localspiketimes(localspikes(j)) * 1000;
        startx = st+eit(1);
        width = eit(end) - eit(1);
        recth = rectangle('Position', [startx ylim(1) width diff(ylim)], ...
            'FaceColor', 0.85*[1 1 1], 'EdgeColor', 'none');
    end
    hold on;
    
    trace = rawdata(:, elecs(i)+1);
    trace = trace - mean(trace);
    plot(rawt, trace, 'k');
    axis tight;
    set(gca, 'yLim', [-300 250]);
    ylabel(sprintf('%d', i));

end
set(rawax(1:end-1), 'XTickLabel', []);
xlabel ms

% Plot EI traces
ylims = {[-300 250] [-70 40] [-70 40] [-70 40] [-70 40] [-70 40]};
for i = 1:length(elecs)
    eitraceax(i) = sanesubplot(1+length(elecs), 8, {1+i 4});
    trace = ei(elecs(i),eiidx);
    trace = trace - mean(trace);
    plot(eit, trace, 'k');
    axis tight;
    set(eitraceax(i), 'YLim', ylims{i});
end
set(eitraceax(1:end-1), 'XTickLabel', []);
xlabel ms

% Plot EI summary
subplot(1, 2, 2);
plot_electrodes(d.d01.ei.position, 'max_scale', 0.05, 'pos_color', 0.25 * [1 1 1]);
plot_electrodes(d.d01.ei.position(elecs,:), 'max_scale', 0, 'label', true, 'zlevel', 20);
hold on;
plot_ei_(ei, d.d01.ei.position, 0, 'scale', 2, 'zlevel', 10);
axis off;


% For EPS, switch to axis specific rects
% % Highlights
% ylim1 = get(rawax(1), 'YLim');
% figpos1 = dsxy2figxy(rawax(1),   [0 ylim1(2) 0 0]);
% ylim2 = get(rawax(end), 'YLim');
% figpos2 = dsxy2figxy(rawax(end), [0 ylim2(1) 0 0]);
% y = figpos2(2);
% h = figpos1(2) - figpos2(2);


%% EI figure
% Not working because midget spikes too small, swamped by parasols
i = 1;
points = 800;
elecs = [163 175 149 340 407 333];

cid = d.d01.offM(i);
cnum = get_cell_indices(d.d01, cid);
spikes = d.d01.spikes{cnum};


rawpath = fullfile(server_data_path, 'rawdata', piece, 'data001.bin');
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(rawpath);
startspike = 13;
start = (spikes(startspike) - 2e-3) * d.d01.sampling_rate;
rawdata = rawFile.getData(start, points);
rawFile.close();
t = spikes(startspike) + ((1:points) - 1)/d.d01.sampling_rate;

clf();
set(gcf(), 'Color', 'w');

for i = 1:length(elecs)
    sanesubplot(length(elecs), 2, {i 1});
    trace = rawdata(:, elecs(i)+1);
    trace = trace - mean(trace);
    plot(t, trace, 'k');
    axis tight;
    set(gca, 'yLim', [-150 150]);
end

subplot(1, 2, 2);
plot_electrodes(d.d01.ei.position, 'max_scale', 0.05, 'pos_color', 0.25 * [1 1 1]);
hold on;
plot_ei_(get_ei(d.d01, cid), d.d01.ei.position, 0, 'scale', 2, 'zlevel', 10);
axis off;


%% Finetuned
i = 1;
clf();
subplot(2,1,1);
plot_ei_(get_ei(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 2);
subplot(2,1,2);
plot_ei_(get_csd(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 2);

%%
i = 3;
clf();
subplot(2,1,1);
plot_ei_(get_ei(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 2);
subplot(2,1,2);
plot_ei_(get_csd(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 2);

%%
i = 4;
clf();
subplot(2,1,1);
plot_ei_(get_ei(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 1.5);
subplot(2,1,2);
plot_ei_(get_csd(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 1.5);

%%
i = 6;
clf();
subplot(2,1,1);
plot_ei_(get_ei(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 2);
subplot(2,1,2);
plot_ei_(get_csd(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 2);


%%
i = 8;
clf();
subplot(2,1,1);
plot_ei_(get_ei(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 2);
subplot(2,1,2);
plot_ei_(get_csd(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 2);


%%
i = 10;
clf();
subplot(2,1,1);
plot_ei_(get_ei(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 1.5);
subplot(2,1,2);
plot_ei_(get_csd(d.d01, d.d01.offM(i)), d.d01.ei.position, 0, 'scale', 2.5);


%% Parasols
i = 3;
clf();
subplot(2,1,1);
plot_ei_(get_ei(d.d01, d.d01.offP(i)), d.d01.ei.position, 0, 'scale', 2);
subplot(2,1,2);
plot_ei_(get_csd(d.d01, d.d01.offP(i)), d.d01.ei.position, 0, 'scale', 1.5);



%% Streamed runs and cell picks
d.d01s = load_data([piece '/data001'], staopts);
d.d01s.onP = [336 3923 5509 6651];
d.d01s.onM = [2617 2977];

d.d03s = load_data([piece '/data003'], staopts);
d.d03s.onP = [3454 3739 5209 6049 6499 7237];
d.d03s.onM = [367 1011 2614 2855 3064 3919 4597 4971 5543 7042 7127];
d.d03s.onMsearched = 3650;
d.d03s.onMsurround = 3933;
d.d03s.offPsearched = [2521 2971 3616 722 4801 1459 887 2941 4456 3841 3647 4172];
d.d03s.offMsearched = [3152, 2851, 7231, 6931, 856, 5206, 5071, 3646, 3631, 346, 77, 2596, 4591, 4336, 3677];
d.d03s.offMsurround = [3031, 2866, 2776, 2686, 2492, 1336, 1186, 7081, 7351, 7576, 7561, 6706, 3031, 2492, 1186, 452, 451, 332, 5732, 5536, 5522, 5326, 5266, 5056, 4846, 4742, 4576, 3811, 3826, 3721, 3676, 3226, 7698, 7681, 7547 7576, 7561, 138, 122, 48, 827, 556, 541, 1111, 1051, 466, 751, 7576, 1382, 1291, 1217, 961, 1741, 1546, 1336, 1276, 1186, 827, 556, 541, 452, 451, 332, 1488, 1353, 1231, 1081, 3226, 4742, 4186, 5056, 4846, 4742, 4576, 4516 4339, 4337, 4306, 4051, 4066, 4576, 4516, 4411, 4339, 4337, 4306, 4051, 3826, 3721, 3691, 3676];


%        OFF parasols surround: 2102 & 2416 (not isolated)
%        OFF parasol surround: 3902 & 3841 & 3647(have to check cone-by-cone)
%        OFF parasols surround: 1501 (well isolated, not strong)
%        OFF parasols surround: 7726 (well isolated but not that strong), 691 & 7277 & 4188 (not isolated)
%        OFF parasol surround: 4188 fairly well isolated
%        OFF parasol surround: 3902, also surround from ON parasol 3739


%% Interdigitation
d.d01s = get_sta_summaries(d.d01s, {4}, 'robust_std_method', 5);
d.d03s = get_sta_summaries(d.d03s, {4}, 'robust_std_method', 5);
plot_rf_coloring(d.d01s, {4});
plot_rf_coloring(d.d03s, {4});


%% Likelihood surface figure
d.d01s = load_cones(d.d01s, 'bayes');
d.d03s = load_cones(d.d03s, 'bayes');

% Combining blurred things a little too much for my taste; just use d01
% Combine?
% expnt = 0.25;
% subplot(131);
% imagesc(matrix_scaled_up((norm_image(d.d01s.cones.dll)-0.5).^expnt,3));
% axis equal square
% subplot(132);
% imagesc(matrix_scaled_up((norm_image(d.d03s.cones.dll)-0.5).^expnt,3));
% axis equal square
% subplot(133);
% imagesc(matrix_scaled_up((norm_image(d.d01s.cones.dll + d.d03s.cones.dll)-0.5).^expnt,3));
% axis equal square

expnt = 0.25;
imagesc(matrix_scaled_up((norm_image(d.d01s.cones.dll)-0.5).^expnt,3));
axis equal square

% 100 um scale bar, based on Greg's measurement of the rig A OLED total screen size
plot([2500 2500-8.8235*3*10], [2600 2600], 'w')


%% data006cm This is actually good for searching too...
% Didn't include ON midget 2614 in original search; should go back and look
% at parasols for this one.
d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], loadopts);

d.d06cm = read_stim_lisp_output(d.d06cm);
d.d06cm.stimulus = parse_stim_rgbs(d.d06cm.stimulus);
urgb = d.d06cm.stimulus.urgb;

cell_list_map_d03s = map_ei(d.d03s, d.d03_06_07);

d.d03s = load_sta(d.d03s, 'load_sta', []);
d.d03s = set_polarities(d.d03s);
d.d03s = get_sta_fits_from_vision(d.d03s);
singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

% Search all ON Parasols
celltype = 1;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = d.d03s.cell_types{celltype}.cell_ids
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false);
end
% Interesting one found!     2988 (matching ON midget: 2855, matching OFF parasols: )
% Slightly interesting ones: 2768 3920

% Search all OFF Parasols
celltype = 2;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = d.d03s.cell_types{celltype}.cell_ids
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Lots of interesting, if hard to work with, stuff in here!!
%
% Lots of costimulated, good number of surround stimulated too!  But not
% many good cases of isolated stimulation, especially surround stimulation
%
% Especially nice ones, with overlap on ON midgets or even ON parasols:
%   ON parasol 2988 (only fires with double cone ON)
%        ON midget 2855 is pretty,
%        OFF parasols: 2521 (reasonably isolated), 2971 (not isolated)
%        OFF parasols surround: 2102 & 2416 (not isolated)
%   ON parasol 3739
%        OFF parasol: 3616
%        OFF parasol surround: 3902 & 3841 & 3647(have to check cone-by-cone)
%   ON parasol 7237
%        OFF parasol: 722 (not completely isolated)
%   ON parasol 5209
%        OFF parasol: 4801 (not completely isolated, check cone-by-cone)
%   ON midget 1011
%        OFF parasols: 1459 (reasonably isolated, edge), 887 (not completely isolated)
%        OFF parasols surround: 1501 (well isolated, not strong)
%   ON midget 367
%        OFF parasols surround: 7726 (well isolated but not that strong), 691 & 7277 & 4188 (not isolated)
%   ON midget 3064
%        OFF parasol edge: 2941 (not well isolated; interesting to look cone-by-cone)
%   ON midget 4597
%        OFF parasol: 4456 well isolated
%        OFF parasol surround: 4188 fairly well isolated
%   ON midget surround 3933
%        OFF parasol: 3841 & 3647, also surround from ON parasol 3739; check cone-by-cone
%        OFF parasol edge: 4172, also surround from ON midget 4597 and others
%        OFF parasol surround: 3902, also surround from ON parasol 3739

% Search all ON Midgets
celltype = 3;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = d.d03s.cell_types{celltype}.cell_ids
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Not too much here; mostly just the intended ones hit
%   ON parasol 3739
%       ON midget 3650 pretty decent for a few of the cones
%   ON midget 6802, 6197, very weak
%   ON midget 3336 not really any good
%
% ON midget 3933 shows surround stimulation so that's cool; should search back for the region in OFF parasols.

% Search all OFF Midgets
celltype = 4;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = d.d03s.cell_types{celltype}.cell_ids(end)
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 15);
end
% Hoo boy...  175 mapped OFF midgets (201 in the streamed data).  Lots of
% crazy stuff in here.  Definitely need to check spike sorting...
%
%   ON parasol 2988 & ON midget 2855
%       OFF midget edge: 3152, 2851, 2641
%       OFF midget surrounds: 2401
%       OFF midget surrounds shared: 3031, 2866, 2776, 2686, 2492, 2341, 2251, 1741, 1546, 1336, 1186 %%
%   ON parasol 7237
%       OFF midget center: 7231, 6931
%       OFF midget edge: 856
%       OFF midget surrounds (may have mixed up some shared in here): 7081, 7351, 7336, 7576, 7561, 6706
%       OFF midget surrounds shared: 3031, 2492, 1186, 452, 451, 332
%   ON parasol 5209
%       OFF midget center: 5206
%       OFF midget edge: 5191, 5071
%       OFF midget surrounds: 5732, 5536, 5522, 5326, 5266
%       OFF midget surrounds shared: 5056, 4846, 4742, 4576
%   ON parasol 3739
%       OFF midget center+surround shared: 3646, 3631
%       OFF midget surrounds: 3811
%       OFF midget surrounds shared: 3826, 3721, 3676, 3226
%   ON midget 367
%       OFF midget center: 346 (perfect for U/D)
%       OFF midget edge: 421, 77
%       OFF midget surrounds: 7698, 7681, 7547 (well isolated), 7576, 7561, 138, 122, 48 %%
%       OFF midget surrounds shared: 827, 556, 541
%       OFF midget, not much despite being closer: 7756
%       OFF midget, not much despite being seemingly close enough: 7516, 751 %%
%   ON midget 1011
%       OFF midget edge: 1111, 1051, 466
%       OFF midget edge shared: 751
%       OFF midget surrounds: 7576, 1382, 1291, 1217, 961
%       OFF midget surrounds shared: 1741, 1546, 1336, 1276, 1186, 827, 556, 541, 452, 451, 332 %%
%       OFF midget weak surround: 1488, 1353, 1231, 1081
%   ON midget 2614
%       OFF midget center: 2596 (maybe usable for OFF midget U/D)
%   ON midget 3064
%       OFF midget edge: 3211
%       OFF midget surrounds shared: 3226
%   ON midget 4597
%       OFF midget center: 4591
%       OFF midget edge: 4336
%       OFF midget surrounds: 4742, 4186
%       OFF midget surrounds shared: 5056, 4846, 4742, 4576, 4516 4339, 4337, 4306, 4051 %%
%   ON midget surround 3933
%       OFF midget center: 3677 (perfect; can use for OFF U/D)
%       OFF midget edge: 4141, 4036
%       OFF midget surround: 4066
%       OFF midget surrounds shared: 4576, 4516, 4411, 4339, 4337, 4306, 4051, 3826, 3721, 3691, 3676 %%
%   Nulls, i.e. perhaps too far for surround stim (but could also just be shitty cell...)
%       OFF midget 1951 (were couple more just after this), 1786, 1771, 1756 %%

% SBCs
celltype = 5;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = d.d03s.cell_types{celltype}.cell_ids
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 15);
end


%% data006cm Revisit search using mapping from d05s; in particular, many more ON midgets there
d.d03_06_07 = load_data([piece '/d03-06-07-norefit/d03-06-07-norefit'], loadopts);
d.d06cm = load_data([piece '/d03-06-07-norefit/data006/data006'], loadopts);
d.d07cm = load_data([piece '/d03-06-07-norefit/data007/data007'], loadopts);

d.d05s = load_data([piece '/data005'], loadopts);
cell_list_map_d05s = map_ei(d.d05s, d.d03_06_07);

d.d06cm = read_stim_lisp_output(d.d06cm);
d.d06cm.stimulus = parse_stim_rgbs(d.d06cm.stimulus);
urgb = d.d06cm.stimulus.urgb;

d.d05s = load_sta(d.d05s, 'load_sta', []);
d.d05s = set_polarities(d.d05s);
d.d05s = get_sta_fits_from_vision(d.d05s);

singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

conerun = d.d05s;
rasterrun = d.d06cm;

% Not much for ON Parasols

% Search all OFF Parasols
celltype = 2;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Seems like some new ones in here:
%   ON parasol 2988 & ON midget 2855 (d03s)
%       OFF parasol shared surround 1907 (d05s)
%
%   ON parasol 3739 (d03s)
%       OFF parasol edge 3616 (d05s)
%       OFF parasol shared edge 3857 (d05s)
%       OFF parasol shared surround weak 3901 (d05s)
%
%   ON parasol 5209 (d03s super weak)
%       OFF parasol shared edge 4666 (d05s)
%
%   ON parasol 7237 (d03s)
%       OFF parasol edge 6721, 7381 (d05s)
%       OFF parasol shared surround 7126 (d05s)
%
%   ON midget 367 (d03s)
%       OFF parasol center 122 (d05s)
%       OFF parasol edge 422 (d05s)
%       OFF parasol shared surround 766, 7126 (d05s)
%       OFF parasol isolated surround 7726 (d05s)
%
%   ON midget 1011 (d03s)
%       OFF parasol edge 842, 1277 (d05s)
%       OFF parasol shared surround 766, 1907 (d05s)
%
%   ON midget 2614 (d03s)
%       OFF parasol shared surround 1907 (d05s)
%
%   ON midget 3933 (d03s)
%       OFF parasol edge 4261 (d05s)
%       OFF parasol shared edge 3857 (d05s)
%       OFF parasol shared surround weak 3901 (d05s)
%
%   ON midget 4597 (d03s)
%       OFF parasol center 4456 (d05s)
%       OFF parasol shared edge 4666 (d05s)
%
%   OFF parasol center 5146, 5941, 6091 (d05s)
%   OFF parasol edge 4352, 6574 (d05s)


% Search all ON Midgets
celltype = 3;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Couple new ones
%   ON parasol 3739 (d03s)
%       ON midget edge 3738 (d05s)
%       ON midget shared surround 3706 (d05s)
%
%   ON midget 3933 (d03s)
%       ON midget shared surround 3706 (d05s)
%
% ON midget edge 6346 (d05s)
% ON midget edge weak 6803 (d05s)
% ON midget edge 3650 (d03s; already found) -> 3648 d05s
% ON midget shared surround 3933 (d03s; already found) -> 3931 d05s
%


%% Condensed rasters for specific interesting groups of cells
% (pick set here and then run next cell to plot)

% All four types (weak ON parasol)
ons = [2988 2855];
ons_d05s = [];
offs = [2521 2416 2102 2851 2341];
offs_d05s = [];
%offs = [2521 2971 2102 2416]; % P
%offs = [offs 3152, 2851, 2641, 2401]; % M
%offs = [offs 3031, 2866, 2776, 2686, 2492, 2341, 2251, 1741, 1546, 1336, 1186]; % M shared

% All four types (weak ON parasol)
ons = [3739 3650];
ons_d05s = [3738 3706];
offs = [3616 3631     3646 3902 3826 3811];
offs_d05s = [];
%offs = [3616 3902 3841 3647   3646 3631 3811 3826 3721 3676 3226];

% Divergence
% ons = []; ons_d05s = []; offs = []; offs_d05s = [];
ons = [4597]; ons_d05s = []; offs = [4591]; offs_d05s = []; % ON and OFF M, only one strong cone in OFF, but looks pretty good.
ons = []; ons_d05s = []; offs = [3677 4172]; offs_d05s = []; % OFF M and P, pretty decent
% ons = [367]; ons_d05s = []; offs = [346 77]; offs_d05s = []; % ON and OFF M, but ON M only ~3 Hz.
% ons = [1011]; ons_d05s = []; offs = [1111 1051 466 751]; offs_d05s = [842 1277];
% ons = [2614]; ons_d05s = []; offs = [2596]; offs_d05s = [1907];
% ons = [3064]; ons_d05s = []; offs = [3211 2941]; offs_d05s = [];


% Surprising surround effects
ons = [367];
ons_d05s = [];
offs = [7726 691 7277          346 77         48 421 122 138   7756   7561 7681 7698 7547 7516];
offs_d05s = [];
%offs = [7726 691 7277 4188            346 421, 77 7698, 7681, 7547, 7576, 7561, 138, 122, 48  827, 556, 541  7756 7516 751];

% Nice surround gradient
ons = [1011];
ons_d05s = [];
offs = [1051   1111, 466    1081, 1276, 1217, 961, 1353, 1291, 1382, 1546, 1741, 1231, 1488];
offs_d05s = [];
%offs = [1111, 1051, 466       751         7576, 1382, 1291, 1217, 961       1741, 1546, 1336, 1276, 1186, 827, 556, 541, 452, 451, 332     1488, 1353, 1231, 1081];
% Look back through parasols

% ON midget surround
ons = [3933];
ons_d05s = [];
offs = []; % FILL IN
offs_d05s = [];


%% Condensed rasters for specific interesting groups of cells
% (pick set above and then run this to plot)
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = ons;
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2, 'hist_line_width', 2, 'hist_color', 'm');
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnycpoly{1}, 'fit', false, 'az_pad_factor', 10);
    set(gcf, 'Position', [516         186        1188         420]);
end
for cellid = ons_d05s;
    cellnum = d.d05s.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2, 'hist_line_width', 2, 'hist_color', 'm');
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d05s, cellid, d.d06cm.stimulus.mapnycpoly{1}, 'fit', false, 'az_pad_factor', 10);
    set(gcf, 'Position', [516         186        1188         420]);
end
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = offs;
    cellnum = d.d03s.cell_nums(cellid);
    mappedcellid = cell_list_map_d03s{cellnum};
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2, 'hist_line_width', 2, 'hist_color', 'm');
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d03s, cellid, d.d06cm.stimulus.mapnycpoly{1}, 'fit', false, 'az_pad_factor', 10);
    set(gcf, 'Position', [516         186        1188         420]);
end
for cellid = offs_d05s;
    cellnum = d.d05s.cell_nums(cellid);
    mappedcellid = cell_list_map_d05s{cellnum};
    figure
    urgb_raster_subplot(d.d06cm, mappedcellid, d.d06cm.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2, 'hist_line_width', 2, 'hist_color', 'm');
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(d.d05s, cellid, d.d06cm.stimulus.mapnycpoly{1}, 'fit', false, 'az_pad_factor', 10);
    set(gcf, 'Position', [516         186        1188         420]);
end


%% data004cm; do some searching here too
d.d01_04_05 = load_data([piece '/d01-04-05-norefit/d01-04-05-norefit'], loadopts);
d.d04cm = load_data([piece '/d01-04-05-norefit/data004/data004'], loadopts);
d.d05cm = load_data([piece '/d01-04-05-norefit/data005/data005'], loadopts);

cell_list_map_d01s = map_ei(d.d01s, d.d01_04_05);

d.d04cm = read_stim_lisp_output(d.d04cm);
d.d04cm.stimulus = parse_stim_rgbs(d.d04cm.stimulus);
urgb = d.d04cm.stimulus.urgb;

d.d01s = load_sta(d.d01s, 'load_sta', []);
d.d01s = set_polarities(d.d01s);
d.d01s = get_sta_fits_from_vision(d.d01s);
singleon  = find(urgb.singles & any(urgb.incr) & any(urgb.absolute_intensities == 1.44));
singleoff = find(urgb.singles & any(urgb.decr) & any(urgb.absolute_intensities == 1.44));
updown    = find(urgb.doubles & any(urgb.incr) & any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleon  = find(urgb.doubles & ~any(urgb.decr) & sum(urgb.absolute_intensities) == 2.88);
doubleoff = find(urgb.doubles & ~any(urgb.incr) & sum(urgb.absolute_intensities) == 2.88);

conerun = d.d01s;
rasterrun = d.d04cm;


% Search all ON Parasols
celltype = 1;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d01s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false);
end
% ON parasol 336 responds to double cones


% Search all OFF Parasols
celltype = 2;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = [7098, 6106, 4246, 3902]%conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d01s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Multiple region OFF parasol surrounds
%   7278, 5990, 5853, 5747, 4801, 3618, 2822, 2748, 2042, 1324
%
% Cusp OFF parasol surround
%   1501, 1381, 46 (ON parasol 336 region)
%
% One region dominant OFF parasol surround
%   7098, 6106, 4246, 3902
%
% Edge OFF parasol
%   5043, 4172, 3724, 2268, 1188 & 261 (ON parasol 336 region)
%
% Mostly isolated OFF parasol center
%   3826, 2521, 887 (inside cusp?)
%
% Non-isolated OFF parasol center
%   6886, 6858, 6423
%
% Not much stimulated OFF parasol
%   6797, 5896, 4591, 692, 151


% Search all ON Midgets
celltype = 3;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d01s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% Edge: 3933
% Center: 2977
% Should search for more


% Search all OFF Midgets
celltype = 4;
ploturgbs = {{singleoff} {updown} {doubleoff}; {singleon} {} {doubleon}};
titles = {'1 cone down' 'down/up' '2 cones down'; '1 cone up' '' '2 cones up'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d01s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 15);
end
% Center OFF midgets
%   466 (ON parasol 336 region)
%   6558, 3781, 3287, 2852, 2596
%
% Cusp OFF midgets
%   2584
%
% Edge OFF midgets
%   335, 332 (ON parasol 336 region)
%   6841, 6664, 6631, 6601, 6511, 6496, 6392, 6391, 6256, 5626, 3676, 2642, 2641, 2431, 2192, 871
%
% Isolated surround OFF midgets:
%   7711, 828, 529, 347, 138, 122 (ON parasol 336 region)
%   707, 423 (ON parasol 336 region weak)
%   6466, 4173, 4127, 4112, 4051, 3056, 3828, 3811, 2776, 2747
%   3797, 3736 (very weak)
%
% Mostly isolated surround OFF midgets:
%   7681 (ON parasol 336 region)
%
% Dominated surround OFF midgets
%   7578, 7561, 7366, 1111, 1081, 1051, 916, 886, 796, 451 (ON parasol 336 region)
%   7396, 7264, 7248, 7247, 7096, 7081, 6991, 6856, 6856, 6811, 6707, 6560, 6272, 5986, 5716, 5582, 5506, 5492, 5237, 5209, 4771, 4411, 4338, 4143, 4142, 3721, 3646, 3632, 3316, 2356, 1741, 1727, 1338 %
%   6946, 6093 (low baseline)
%   5297, 4592, 4186, 1982 (weak)
%   5267, 5073, 5071, 4726, 4516, 4340, 4336, 2028, 1892, 1681, 1231i, 1173 (very weak)
%
% Multiple surround OFF midgets
%   7548, 7471, 7187, 7173, 7158, 7141, 6931, 6917, 6587, 6303, 6302, 6287, 6286, 5611, 4966, 4876, 4711, 4501, 4486, 4233, 4216, 4007, 3692, 3466, 3455, 3406, 3151, 3121, 3106, 3046, 3031, 2866, 2686, 2566, 2492, 2401, 2252, 2191, 1471, 856 %
%   7516, 7501, 6031, 5943, 5881, 5807, 5672, 5416, 4847, 4846, 3616, 3362, 3226, 2386, 1801, 1666, 1531, 1516, 1277 (weak) (??? possibly mislabeled)
%
% Not much stimulated
%   5671, 5657, 5401, 5329, 5327, 5191, 5056, 4398, 3061, 1891, 1771, 1756, 1637, 1486, 1426, 1367
%   2056, 1307 (fairly close)
%   961, 49, 34 (fairly close ON parasol 336 region)
%
% Weird OFF midgets
%   7232, 6062, 5522, 5178


% Search all SBCs
celltype = 5;
ploturgbs = {{singleon} {updown} {doubleon}; {singleoff} {} {doubleoff}};
titles = {'1 cone up' 'up/down' '2 cones up'; '1 cone down' '' '2 cones down'};
for cellid = conerun.cell_types{celltype}.cell_ids
    cellnum = conerun.cell_nums(cellid);
    mappedcellid = cell_list_map_d01s{cellnum};
    if isempty(mappedcellid), continue; end
    
    figure
    urgb_raster_subplot(rasterrun, mappedcellid, rasterrun.triggers(1:2:end), ploturgbs, 'titles', titles, 'add_width', 2);
    
    sanesubplot(2, 5, {1:2, 4:5});
    plot_rf_stimmap(conerun, cellid, rasterrun.stimulus.mapnyc{1}, 'fit', false, 'az_pad_factor', 10);
end
% 2313 shows some nice surround.