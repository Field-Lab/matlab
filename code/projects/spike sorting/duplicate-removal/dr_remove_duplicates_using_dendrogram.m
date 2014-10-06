% dr_remove_duplicates_using_dendrogram
%
% this is a complete duplicate finding algorithm that uses a distance
% function (either spike times or eis are used for comparisons between
% cells) and then it constructs a single linkage dendrogram and thresholds
% it. each duplicate group found in the dendrogram is merged (ie no cells
% are discarded).
%
% this is probably the best starting point for looking at duplicate finding
% algorithms/distance metrics.
%
tic

clear all

% set up parameters
THRESHOLD      = .65;     % threshold for choosing duplicate classes
SAMPLING_RATE  = 20000;   % sampling rate (20 khz)
N_ELECTRODES   = 61;      % how many electrodes in the dataset
USE_EIS        = true;    % use spike times or eis
UPDATE_PATH    = false;   % read from a copy of the neurons file or not
SPACE_ONLY_EIS = false;   % use spatial or spatiotemporal eis

% begin
disp(sprintf('Duplicate Removal: Threshold set to %0.2f',THRESHOLD))

% initialize datarun
% caution: when duplicate removal has finished, the neurons file at this
% path will be updated with the cleaned neurons (ie, the raw neurons file
% will be changed)
datarun = load_data('/Data/Machado/2009-04-13-0/data000/test-dr/test-dr');

% update the path to the neurons-raw file:
if UPDATE_PATH
    % name of the neurons file (that will be modified)
    datarun.names.rrs_neurons_path = [datarun.names.rrs_neurons_path '-copy'];
    % name of input ei file (computed on ALL cells in neurons-raw file)
    datarun.names.rrs_ei_path = [datarun.names.rrs_ei_path];
end

% load neurons
datarun = load_neurons(datarun);
indices = get_cell_indices(datarun, datarun.cell_ids);

% compare spike trains or eis
if USE_EIS
    % load eis
    datarun = load_ei(datarun, 'all', 'array_type', N_ELECTRODES);
    % remove duplicates
    [corr_result] = dr_ei_duplicate_removal(datarun,'verbose',1,'corr_threshold',THRESHOLD,'space_only',SPACE_ONLY_EIS);

else
    % compute distance matrix
    disp('Computing Distance Matrix...')
    spiketrain=zeros(datarun.duration*1000/10,length(datarun.spikes)+length(datarun.spikes));
    for i=1:length(datarun.spikes)
        if isempty(datarun.spikes{i}), continue; end
        spiketrain(:,i)=histc(datarun.spikes{i},0:10/1000:datarun.duration-10/1000);
    end
    for i=1:length(datarun.spikes)
        if isempty(datarun.spikes{i}), continue; end
        spiketrain(:,i+length(datarun.spikes))=histc(datarun.spikes{i},0:10/1000:datarun.duration-10/1000);
    end

    corr_result = corrcoef(spiketrain);
    corr_result = corr_result(1:length(datarun.spikes),length(datarun.spikes)+1:end);
    corr_result = setdiag(corr_result,1);
end

corr_result = 1-abs(corr_result);
THRESHOLD   = 1-THRESHOLD;

% perform complete-linkage clustering since we only want to merge maximal cliques 
disp('Constructing Linkages...')
corr_links = linkage(squareform(corr_result), 'complete');

% cluster at specified threshold
clusters = cluster(corr_links, 'Cutoff', THRESHOLD, 'Criterion', 'distance');

% plot it
figure; dendrogram(corr_links, 0, 'colorthreshold', THRESHOLD);
figure; hist(corr_links(:,3),100)

% merge all spike trains in each cluster
unique_cells.s = {};
unique_cells.e = [];

for cluster = 1:length(datarun.spikes)
    datarun.spikes{cluster} = datarun.spikes{cluster} * SAMPLING_RATE;
end

disp('Merging Spike Trains...')
for cluster = 1:length(unique(clusters))

    ids = find(clusters == cluster);
    merged_cell.s = datarun.spikes{ids(1)};
    merged_cell.e = datarun.channels(ids(1));
    
    if length(ids) > 1
        for current_id = 2:length(ids)
            current_cell.s = datarun.spikes{ids(current_id)};
            current_cell.e = datarun.channels(ids(current_id));
            s.a = current_cell.s; s.b = merged_cell.s;
            e.a = current_cell.e; e.b = merged_cell.e;
            if isempty(s.a) || isempty(s.b), continue; end
            merged_cell = dr_merge(s, e);
        end
    end
    
    unique_cells.s{end + 1} = merged_cell.s; %#ok<AGROW>
    unique_cells.e(end + 1) = merged_cell.e; %#ok<AGROW>
end

% save out a neurons file
disp(sprintf('%d Unique Cells Found!', length(unique_cells.s)));

% open neurons file
nf = edu.ucsc.neurobiology.vision.io.NeuronFile(datarun.names.rrs_neurons_path);

% get rid of all neurons
kill = [];
for ii = 1:length(datarun.cell_ids)
    nf.deleteNeuron(datarun.cell_ids(ii));
end

id = 0;
for ii = 1:length(unique_cells.s)
    % add nnneurons
    il = edu.ucsc.neurobiology.vision.util.IntegerList;
    spikes = unique_cells.s{ii};
    for jj = 1:length(unique_cells.s{ii})
        il.add(spikes(jj));
    end

    nf.addNeuron(unique_cells.e(ii), id, il,length(unique_cells.s{ii}))
    id = id + 1;
end

nf.close;
toc