function [sta_fits] = load_snls(server_path, cell_ids, channels)
% LOAD_SNL     Loads SNL from Obvius files
%
% usage:   snl = load_snl(server_path, cell_ids, channels);
%
% arguments:   server_path -- path to STA fits directory
%                 cell_ids -- Nx1 vector of *all* cell ID's
%                 channels -- Nx1 vector of *all* electrodes
%
% outputs:        snl -- Nx1 cell array of STA fits
%
% Loads snl from Obvius files. Used internally by load_series.
%
% snl = load_snls('/snle/lab/Experiments/Array/Analysis/2005-08-08-0/rf-10-mg-4.5/', extras.cell_ids, extras.channels);
%
%
% greschner 2008-02-01 from LOAD_STA_FITS
%
    

% determine sta direcotry
sta_dir = dir([server_path '/stas-*']);

% find sta directory
if (length(sta_dir) ~= 1) | (sta_dir(1).isdir ~= 1)
  error('unable to find sta directory');
end

% create new server path
server_path = [server_path '/' sta_dir.name '/'];
files = dir([server_path 'snl-*']);

% allocate memory
sta_fits = cell(length(files),1);
sta_ids  = zeros(length(files),1);

% loop through files
for j=1:length(files)
  
  % check if file exists 
  if ~exist([server_path files(j).name], 'file');
    error(['load_index: can not open ' full_path files(j).name]);
  end

  % load sta fit
  sta_fit = load_snl(server_path, files(j).name);

  % determine electrode
  electrode = find_electrode(files(j).name);
  letter = files(j).name(end);
  
  % find neurons associated with electrode
  indices = find(channels == electrode);
  
  % find cell ID's associated with channel
  candidate_IDs = cell_ids(indices);
  
  % determine specific cell ID associated with electrode
  switch letter
   case 'A', cell_id = candidate_IDs(1);
   case 'B', cell_id = candidate_IDs(2);
   case 'C', cell_id = candidate_IDs(3);
   case 'D', cell_id = candidate_IDs(4);
   case 'E', cell_id = candidate_IDs(5);
   case 'F', cell_id = candidate_IDs(6);
   otherwise, error(['bad letter identifier: ' files(j).name]);
  end
  
  % save cell ID number
  sta_fit.cell_id = cell_id;

  % save all sta fits
  sta_fits{j} = sta_fit;
  
end

% this is good insurance and it deals with neurons which have no
% sta fits by assigning blanks
sta_fits = sort_by_ids(sta_fits, cell_ids);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sort_by_ids(in, ids)
% sort STA fits by the cell IDs

out = cell(size(in));
in_IDs = [];

% check for no stas
no_sta = length(ids) - length(in);
if (no_sta > 0)
  disp(['load_snl: no snl for ' num2str(no_sta) ' cells.']);
end

% extra "in" ID's
for i=1:length(in)
  in_IDs = [in_IDs in{i}.cell_id];
end

% re-order items in cell array
for i=1:length(ids)
  index = find(ids(i) == in_IDs);
  % check if fit exists
  if ~isempty(index)
    % re-order correctly
    out{i} = in{index};
  else
    % otherwise, assign blank
    out{i} = [];
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sta_fit = load_snl(server_path, file_name)
% load individual Obvius STA fit

% load file data
fid = fopen([server_path file_name]); 

% read in all text
file=textscan(fid, '%s'); file_text=file{1};
fclose(fid);

% loop through entries in text
i=0;
while i<length(file_text)
  i=i+1;
  
  % found new field name
  if ~isempty(strfind(file_text{i},':'));
    
    % format identified field name 
    field_name = strrep(file_text{i},':','');
    field_name = strrep(field_name,'(','');
    field_name = strrep(field_name,')','');
    field_name = strrep(field_name,'-','_');
    field_name = lower(field_name);

    % create blank structure field name
    sta_fit.(field_name)=[];
    
  else
  
    % format identified field value
    field_value = strrep(file_text{i},'(','');
    field_value = strrep(field_value,')','');
    field_value = strrep(field_value,'#','');

    % set identified field value
    if ~isempty(str2num(field_value))
      sta_fit.(field_name)=[sta_fit.(field_name) str2num(field_value)];
    else
      sta_fit.(field_name)=[sta_fit.(field_name) field_value];
    end
  
  end
end 

% set file name
sta_fit.obvius_id=file_name;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function electrode = find_electrode(file_name)
% determine the electrode name

% if 512 electrode array
if file_name(5) == 'L'
  electrode = str2num(file_name(6:8));
  return
end

% if 61 electrode array
%error('61 not tested')

electrode_name = file_name(5:6);
electrode = convert_electrode_name(electrode_name);







