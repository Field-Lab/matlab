function dat = loadData(computer,dataSet)

if ~exist('computer','var') || isempty(computer)
  computer = 'stanford';
end

switch computer
  case 'lcv'
    mainDataPath = '/share/users-lcv/freeman/retinalData/data/';
    mainAnalPath = '/share//users-lcv/freeman/retinalData/anal/';
    mainFigurePath = '/share/users-lcv/freeman/retinalData/figure/';
  case 'laptop'
    mainPath = '~/Dropbox/Projects/retinaSubunits/';
    mainDataPath = '~/Dropbox/Projects/retinaSubunits/data/';
    mainAnalPath = '~/Dropbox/Projects/retinaSubunits/anal/';
    mainFigurePath = '~/Dropbox/Projects/retinaSubunits/figure/';
  case 'lcv-2'
    mainDataPath = '/share/users-lcv/freeman/retinalData/data-2/';
    mainAnalPath = '/share//users-lcv/freeman/retinalData/anal-2/';
    mainFigurePath = '/share/users-lcv/freeman/retinalData/figure-2/';
  case 'laptop-2'
    mainPath = '~/Dropbox/Projects/retinaSubunits/';
    mainDataPath = '~/Dropbox/Projects/retinaSubunits/data-2/';
    mainAnalPath = '~/Dropbox/Projects/retinaSubunits/anal-2/';
    mainFigurePath = '~/Dropbox/Projects/retinaSubunits/figure-2/';
  case 'salk'
    mainDataPath = '/snle/acquisition/freeman/';
    mainAnalPath = '/Users/vision/Desktop/freeman/anal/';
    mainFigurePath = '/Users/vision/Desktop/freeman/figure/';
  case 'stanford'
    mainPath = '/Users/alexth/Desktop/Freeman/';
    mainDataPath = '/Users/alexth/Desktop/Freeman/data/';
    mainAnalPath = '/Users/alexth/Desktop/Freeman/anal/';
    mainFigurePath = '/Users/alexth/Desktop/Freeman/figure/';
  case 'bertha'
      tmp = regexp(dataSet,'/');
      mainPath = ['/Volumes/Analysis/',dataSet(1:tmp),'subunits',dataSet(tmp:end),'/'];
      mainDataPath = ['/Volumes/Analysis/',dataSet(1:tmp),'subunits',dataSet(tmp:end),'/'];
      mainAnalPath = ['/Volumes/Analysis/',dataSet(1:tmp),'subunits',dataSet(tmp:end),'/anal/'];
      mainFigurePath = ['/Volumes/Analysis/',dataSet(1:tmp),'subunits',dataSet(tmp:end),'/figure/'];
  otherwise
    dat = [];
    fprintf('(loadData) invalid location specified \n');
    return
end

dataPath = mainDataPath;
analPath = mainAnalPath;
figurePath = mainFigurePath;
% dataPath = strcat(mainDataPath,dataSet);
% analPath = strcat(mainAnalPath,dataSet);
% figurePath = strcat(mainFigurePath,dataSet);

if ~isdir(analPath); mkdir(analPath); end
if ~isdir(figurePath); mkdir(figurePath); end
dat.dataSet = dataSet;
dat.mainPath = mainPath;
dat.dataPath = dataPath;
dat.analPath = analPath;
dat.figurePath = figurePath;
dat.computer = computer;
dat.mainDataPath = mainDataPath;
dat.mainAnalPath = mainAnalPath;
dat.mainFigurePath = mainFigurePath;

switch dataSet
  case 'vornoi-gaussian-08-31-2011'
    %generator signal for each cone (1870 cones in all)
    switch computer
      case 'lcv'
        load(fullfile(dataPath,'cone_inputs'));
        chunkSize = size(cone_inputs,1);
        dat.coneInputs = cone_inputs(1:chunkSize,:);
        clear cone_inputs
      case 'laptop'
        load(fullfile(dataPath,'cone_inputs'));
        chunkSize = size(cone_inputs,1);
        dat.coneInputs = cone_inputs;
        clear cone_inputs
    end
    
    %the number of spikes in each time bin (a bin = 1 monitor frame)
    load(fullfile(dataPath,'spike_rate'));
    dat.spikeRate = spike_rate(:,1:chunkSize);
    clear spike_rate
    
    % identifies the type of each RGC (i.e. ON-midget cell)
    load(fullfile(dataPath,'cell_types'));
    dat.cellTypes = cell_types;
    clear cell_types
    
    %%
    % x-y coordinates for each cone.
    load(fullfile(dataPath,'cone_locations'));
    
    [dat.locations scaleVals] = rescaleLocations(cone_locations);
    
    clear cone_locations
    
    %the weight from each cone to each RGC
    load(fullfile(dataPath,'cone_weights'));
    dat.coneWeights = cone_weights;
    clear cone_weights
    
    %numbers identify each RGC
    load(fullfile(dataPath,'rgc_ids'));
    dat.rgcIds = rgc_ids;
    clear rgc_ids
    
    %x-y coordinates for the center of mass of each RGC RF
    load(fullfile(dataPath,'rgc_locations'));
    
    for ii=1:length(rgc_locations)
      if ~isempty(rgc_locations{ii})
        tmp = rgc_locations{ii};
        tmp(1) = tmp(1) - scaleVals.mnx;
        tmp(2) = tmp(2) - scaleVals.mny;
        tmp = tmp/scaleVals.mx;
        tmp = tmp*scaleVals.s;
        dat.rgcLocations{ii} = tmp;
      else
        dat.rgcLocations{ii} = [];
      end
    end
    %
    clear rgc_locations
    
    %identifies the cones as L,M, or S.
    load(fullfile(dataPath,'cone_types'));
    dat.coneTypes = cone_types;
    clear cone_types
    
  case {'vornoi-binary-10-26-2011','vornoi-ternary-10-26-2011'...
      'vornoi-binary-12-13-2011','vornoi-ternary-12-13-2011'}
    
    load(fullfile(dataPath,'datarun'));
    
    % generator signal for each cone
    dat.coneInputs = datarun.movie';
    
    % the number of spikes in each time bin
    dat.spikeRate = datarun.spike_rate;
    
    % get vornoi regions for each cone
    dat.coneMap = datarun.cones.cone_map;
    dat.vornoiMap = datarun.cones.vornoi_map;
    nCones = max(dat.coneMap(:));
    dat.locations = zeros(nCones,2);
    for ii=1:nCones
      [rows cols] = find(dat.coneMap==ii);
      if ~isempty(rows)
        dat.locations(ii,1) = mean(rows);
        dat.locations(ii,2) = mean(cols);
      else
        dat.locations(ii,:) = [];
      end
    end
    
    dat.locations = rescaleLocations(dat.locations);
    
    %numbers identify each RGC
    dat.rgcIds = datarun.cell_ids;
    
    dat.stas = datarun.stas;
    
  otherwise
    
    warning('off');
    load(fullfile(dataPath,'conepreprocess.mat'));
    warning('on');
    
    if exist('outrun','var');
      datarun = outrun; clear outrun;
    end
    
    if strcmp(datarun.cell_types{4}.name,'OFF Midget nc8')
      datarun.cell_types{4}.name = 'OFF Midget';
      datarun.cell_types{9}.name = 'x';
    end 
    
    dat.spikeRate = double(datarun.spike_rate);
    dat.rgcIds = datarun.cell_ids;
    dat.stas.time_courses = datarun.stas.time_courses;
    dat.stas.polarities = datarun.stas.polarities;
    dat.cellTypes = datarun.cell_types;
    dat.rgcLocations = datarun.stas.rf_coms;
    dat.rfs = datarun.stas.rfs;
    
    % for, vornoi experiments stimulus will be a movie
    % (should be nCones x nTimePoints)
    if isfield(datarun,'cone_inputs')
      dat.coneInputs = double(datarun.cone_inputs);
    elseif isfield(datarun,'movie')
      dat.coneInputs = double(datarun.movie');
    end
    
    % for vornoi experiments, there will be a cone map
    % instead of  cone finding
    if isfield(datarun,'cones')
      if isfield(datarun.cones,'cone_map')
        dat.coneMap = datarun.cones.cone_map;
        nCones = max(dat.coneMap(:));
        dat.locations = zeros(nCones,2);
        for ii=1:nCones
          [rows cols] = find(dat.coneMap==ii);
          if ~isempty(rows)
            dat.locations(ii,1) = mean(cols);
            dat.locations(ii,2) = mean(rows);
          else
            dat.locations(ii,:) = [];
          end
        end
        % should rescale... vornoi map is always 600x600,
        % pixelated stimulus was either 200x200 or 300x300
        dat.locations = dat.locations/3;
      else
        dat.coneWeights = datarun.cones.weights;
        dat.coneTypes = datarun.cones.types;
        dat.locations = datarun.cones.centers;
      end
    else % no cones at all, pretend pixels are cones
      [X Y] = meshgrid([1:datarun.stimulus.field_width],[1:datarun.stimulus.field_height]);
      dat.locations = [X(:) Y(:)];
      dat.coneWeights = zeros(size(dat.locations,1),size(dat.rfs,1));
      for ic=1:size(dat.rfs,1)
        dat.coneWeights(:,ic) = dat.rfs{ic}(:);
      end
      dat.coneTypes = repmat('U',size(dat.locations,1),1);
    end
    
    
end

dat.nCones = size(dat.locations,1);

