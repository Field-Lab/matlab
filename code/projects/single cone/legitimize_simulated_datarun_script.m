% This script produces a datarun struct with single cone RFs, where the cone mosaic and RGC sampling
% are specified in an incomplete datarun struct lacking certain fields.  Typically, the input will be the result
% of a simulation.
%
%
% INPUT FORMAT
%
% The script looks for a variable in memory called "datarun_sim".  If this does not exist,
% it loads the mat file specified in datarun_path, where it expects to find a struct named 'datarun' which
% gets renamed to datarun_sim.
% In either case, the struct must contain these fields:
%
%  datarun_sim.cones.centers            - Nx2 matrix with x,y coordinates of each cone
%                   .types              - N-length char vector of cone types (L,M,S,U for unknown)
%                   .rf_fits.centers    - Mx2 matrix of RF center locations
%                   .weights            - NxM matrix, each column gives the cone weights for one RGC
%
% where N = # cones, M = # RGCs
%
%
% OUTPUT FORMAT
%
% Data from datarun_sim will be used to create in a datarun struct containing all the necessary fields.
% The following fields contain meaningful data:
%
%   datarun.cell_ids                - Mx1 double matrix, numbers 1 through M
%   datarun.cell_types{1}.cell_ids	- Mx1 int32 matrix, numbers 1 through M
%   datarun.cell_types{1}.name      - 'simulated'
%
%   datarun.cones.centers           - Nx2 matrix with x,y coordinates
%                .types             - N-length char vector, 'L','M','S', or 'U' for unknown
%                .rgb               - Nx3 matrix withrgb triplet for each cone
%                .weights           - NxM matrix, each column gives the cone weights for one RGC
%                .rf_fits           - Mx1 cell array of structs, each with fields
%                                       center
%                                       center_radius
%                                       center_scale
%                                       surround_radius
%                                       surround_scale
%                                       fit_radius
%                                       error
%  
%  These fields are also created, but contain fabricated data:
%   likelihoods, likelihoods, types_kmeans, rgb, roi
%
% The new creations are saved in a folder with name "new_name".  They are saved as a datarun struct
% and as text files in the standard format.
% 
%
% 2009-02 gauthier
%





% parameters
datarun_path = [/simulation 'simulations/exclusion/'];
new_name = 'simulations/simulation-balanced-15clumped';



% load simulated datarun
if ~exist('datarun_sim','var')
    % assume it is named 'datarun'
    load(datarun_path)
    % rename to datarun_sim
    datarun_sim = datarun;
end

% clear any previous attempts
clear datarun

% add cell id, type info
num_cells = size(datarun_sim.cones.weights,2);
datarun.cell_ids = 1:num_cells;
datarun.cell_types{1}.cell_ids = int32(datarun.cell_ids);
datarun.cell_types{1}.name = 'simulated';


% enter stimulus info
datarun.stimulus.stixel_width = 1;
datarun.stimulus.stixel_height = 1;
datarun.stimulus.field_height = ceil(max(datarun_sim.cones.rf_fits.center(:,2))*1.2);
datarun.stimulus.field_width = ceil(max(datarun_sim.cones.rf_fits.center(:,1))*1.2);
x_offset = datarun.stimulus.field_width/12;
y_offset = datarun.stimulus.field_height/12;


% enter RGC info
for cc=1:num_cells
    datarun.stas.rf_coms{cc} = datarun_sim.cones.rf_fits.center(cc,:)+[x_offset y_offset];
end


% enter cone info
datarun.cones.weights = datarun_sim.cones.weights;
datarun.cones.types = reshape(char(datarun_sim.cones.types),[],1);
datarun.cones.centers = [datarun_sim.cones.centers(:,1)+x_offset datarun_sim.cones.centers(:,2)+y_offset];

% compute fits to each RF
if ~exist('rf_fits','var')
    datarun = fit_cone_rfs(datarun,{1},'verbose',0,'foa_2d',1);
else
    datarun.cones.rf_fits = rf_fits;
end

% fabricate cone info
num_cones = size(datarun.cones.weights,1);
datarun.cones.likelihoods = zeros(num_cones,1);
datarun.cones.types_em = datarun.cones.types;
datarun.cones.types_kmeans = datarun.cones.types;
datarun.cones.rgb = zeros(num_cones,3);
datarun.cones.roi = ones(num_cones,1);

% save results as text files and a datarun struct
save_dir = [single_cone_path new_name '/'];
mkdir(save_dir)
save([save_dir 'datarun'],'datarun')
export_single_cone_data(datarun,'all',[],save_dir,'verbose',1)
