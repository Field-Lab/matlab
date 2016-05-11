%% Copied out of 2007-09-18-4 match paper analysis.
%% Axon velocities, nicer
piece = '2007-09-18-4';
rig = 'A';
loadopts = struct('load_params', true, 'load_neurons', true);
cell_types = {1, 2, 3, 4, 5};
d02n = load_data([piece '/data002-nwpca'], loadopts);
d02n = load_ei(d02n, cell_types);

%% Manually mark bad eis or bad fits, manually add good eis that were missed.

% Quick way to get original cid indices from current goodcid indexed figure 
% windows (see below).
% simpleprint(find(ismember(d02n.cell_types{j}.cell_ids, goodcids{j}(get(0, 'Children')))))

manualgood = {[], [], [], [], []};
manualbad = {
    [3, 5, 8], [], [163], [], []};
manualbadfits = {  % These could be fixed by fixing eistart/end
  [14, 22, 23, 24, 26, 29, 33, 48, 49], ...
  [13, 17, 22, 25, 45, 47], ...
  [29 30 33 45 53 60 73 81 84 89 101 106 117 119 135 140 157 158 192 195 203 217 221], ...
  [9 10 16 24 31 94 97 104], ...,
  [6, 7, 23]
};

%%
cell_types = [1, 2, 3, 4, 5];
goodcids = cell(size(cell_types));
goodeistarts = cell(size(cell_types));
goodeiends = cell(size(cell_types));
for j = 1:length(cell_types)
  cell_type = cell_types(j);
  threshold = -10;
  cids = d02n.cell_types{cell_type}.cell_ids;
  eistarts = zeros(size(cids));
  eiends = eistarts;
  for i = 1:length(cids)
      cid = cids(i);
      ei = get_ei(d02n, cid);
      [eistarts(i) eiends(i)] = get_axon_ei_times(ei, threshold);
  end

  % Back these up a bit.
  eistarts = eistarts - 2;
  eiends = eiends - 3;

  % Throw out bad ones
  badstart = eistarts < 22 | eistarts > 30;
  badlength = eiends - eistarts < 9;
  good = find(~(badstart | badlength));
  
  % Include manual futzing...
  good = setdiff(good, manualbad{j});
  good = setdiff(good, manualbadfits{j});
  good = [good manualgood{j}];
  
  goodcids{j} = cids(good);
  goodeistarts{j} = eistarts(good);
  goodeiends{j} = eiends(good);
end

fits = cell(size(cell_types));
deltas = cell(size(cell_types));
meanvelocities = cell(size(cell_types));

%%
positions = d02n.ei.position;
A0 = -5000;
sigma0 = 50;

for j = 1:length(cell_types)
  cell_type = cell_types(j);
  cids = goodcids{j};
  eistarts = goodeistarts{j};
  eiends = goodeiends{j};  

  fits{j} = cell(size(cids));
  for i = 1:length(cids)
    fprintf('%d of %d\n', i, length(cids));
    cid = cids(i);
    ei = get_ei(d02n, cid);
    eislice = ei(:, eistarts(i):eiends(i));
    fits{j}{i} = gaussian_fit_ei_negatives(eislice, positions, A0, sigma0);
  end

  deltas{j} = cell(size(fits{j}));
  meanvelocities{j} = zeros(size(fits{j}));
  for i = 1:length(fits{j})
    fit = fits{j}{i};
    deltas{j}{i} = sqrt(sum(diff(fit(2:end, :)).^2, 2));
    meanvelocities{j}(i) = mean(deltas{j}{i}) * 0.02;  % microns / m * 20 kHz
  end
end

%%
j = 3;
cell_type = cell_types(j);
cids = goodcids{j};
for i = 1:(length(cids)/2)
  cid = cids(i);
  figure(i);
  plot_ei(d02n, cid, 'coordinates', 'array', 'neg_color', [0.8 0.8 1], ...
          'pos_color', [1, 0.8, 0.8]);
  plot(fits{j}{i}(2:end,1), fits{j}{i}(2:end,2), 'k-o');
  % title(num2str(meanvelocities{j}(i)));
end

%%
x = 0:0.1:2;
c = {'r', 'r--', 'g', 'g--', 'b'};
hold on
for i = 1:length(cell_types)
  b{i} = histc(meanvelocities{cell_types(i)}, x);
  plot(x, b{i} ./ numel(meanvelocities{i}), c{i});
end
legend({'ON Parasols', 'OFF Parasols', 'ON Midgets', 'OFF Midgets', 'SBCs'});
xlabel 'velocity (m / s)'
ylabel 'frequency'
