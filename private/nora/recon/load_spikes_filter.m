function [spikes, movie, STA, cell_names] = load_spikes_filter(stimtype, exp_nm, cell)

% stimtype is string either 'WN' or 'NSEM'
% exp_nm is string like '2012-08-09-3'
% cell is string of cell type name eg 'ONPar' or array of cell ids eg 841 or [841, 5058]

% spikes come out in (cells, trials) 


spike_directory = '/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/';

% Load and process test movie
[StimulusPars, ~] = StimulusParams(exp_nm, stimtype, 'mapPRJ');
[testmovie0, ~, ~]          = loadmoviematfile(exp_nm , stimtype, '8pix_Identity_8pix','testmovie');
movie             = testmovie0{1}.matrix(:,:,StimulusPars.slv.testframes);
movie = movie - mean(movie(:));
clear testmovie0

% Get all cell names
if ischar(cell)
    matfiles = dir([spike_directory exp_nm '/WN_mapPRJ/STA/STAandROI_' cell '*.mat']);
    n_cell_input = 1;
elseif isnumeric(cell)
    n_cell_input = length(cell);
    for i = 1:n_cell_input
        matfiles(i) =  dir([spike_directory exp_nm '/WN_mapPRJ/STA/STAandROI_*_' num2str(cell(i)) '.mat']);
    end
else
    error('Unrecognizable cell input. Cell must be a string like ONPar or a vector of ids')
end

% Check the number of cells
n_cells = length(matfiles);
if n_cells == 0
    error(['No cell files were found in ' spike_directory])
elseif n_cells < n_cell_input
    warning(['Not all cell files were found in ' spike_directory])
end
clear n_cell_input cell

data_missing = [];

% load up spikes and filter for each cell
for i_cell = 1:n_cells
    load([spike_directory exp_nm '/WN_mapPRJ/STA/' matfiles(i_cell).name])
    STA{i_cell} = flip(STAandROI.STA, 3);
    STA{i_cell} = STA{i_cell} - mean(STA{i_cell}(:));
    clear STAandROI
    
    % Process spikes for glm_execute with proper subroutines
    try
        cell_names{i_cell} = matfiles(i_cell).name(11:(end-4));
        load([spike_directory exp_nm '/' stimtype '_mapPRJ/organizedspikes_' cell_names{i_cell} '.mat']);
        spikes(i_cell, :) = subR_createraster(organizedspikes.block, StimulusPars.slv);
    catch
        data_missing = [data_missing i_cell];
    end
end

cells = 1:n_cells;
for i_cell = 1:n_cells
   if any(isnan(STA{i_cell}(:))) || any(data_missing == i_cell)
       cells = cells(cells ~= i_cell);
   end
end
STA = STA(cells);
spikes = spikes(cells,:);
cell_names = cell_names(cells);

end


function raster_spiketimes = subR_createraster(blockedspikes, TestPars)
% AKHeitman 2014-04-14
% Make a raster which takes into account GLM processing
% blocekdspikes: needs
%   .t_sp_withinblock
%
% TestPars needs
%   .fittest_skipseconds
%   .TestBlocks

rasterblocks = TestPars.TestBlocks;
t_start      = TestPars.fittest_skipseconds;

raster_spiketimes = cell(length(rasterblocks),1);

for i_blk = 1 : length(rasterblocks)
	blknum = rasterblocks(i_blk);
	sptimes = blockedspikes.t_sp_withinblock{blknum} - t_start;
	sptimes = sptimes(find(sptimes > 0 ) );
    % HACK NEEDED FOR 2013-10-10-0 and other long runs
    if isfield(TestPars, 'test_skipENDseconds')
        sptimes = sptimes(find(sptimes < (TestPars.test_skipENDseconds - TestPars.fittest_skipseconds - .1)));
    end
    
    raster_spiketimes{i_blk} = sptimes;
end 

end