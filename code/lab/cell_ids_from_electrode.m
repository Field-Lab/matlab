function cell_ids = cell_ids_from_electrode(electrode)
% cell_ids_from_electrode     which cell ids come from a given electrode
%
% usage:  cell_ids = cell_ids_from_electrode(electrode)
%
% arguments:  electrode - 
%
% outputs:     cell_ids - list of cell ids 
%                           if no argument is taken, the result is printed
%
%
% 2009-06  gauthier
%



max_clust = edu.ucsc.neurobiology.vision.util.VisionParams.maxClusters;

cell_ids = max_clust * double(electrode) + (-14:0)';

if nargout == 0
    fprintf('electrode %d: cell ids %d to %d\n',electrode,min(cell_ids),max(cell_ids))
    clear cell_ids
end
