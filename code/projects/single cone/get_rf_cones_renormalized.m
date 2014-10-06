function datarun = get_rf_cones_renormalized(datarun, cell_spec, varargin)
% get_rf_cones_renormalized     Normalize the weights of each RGC
%
% usage:  datarun = my_function(datarun, arg1, <params>)
%
% arguments:  datarun - datarun struct with datarun.cones struct
%           cell_spec - which cells
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in field X
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           false               show output
%
%
%
% 2008-10 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('type','fit', @(x)any(strcmpi(x,{'fit','max','noise std','center sum','center std'})));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;





% BODY OF THE FUNCTION

cell_indices = get_cell_indices(datarun,cell_spec);

% go through each RGC
for cc = 1:length(cell_indices)
    
    cell_index = cell_indices(cc);
    
    the_weights = datarun.cones.weights(:,cell_index);
    
    % skip if weights are all zero
    if all(the_weights == 0)
        continue
    end
    
    switch params.type
        
        case 'fit' % normalize to the gaussian fit height
            
            % get the fit
            the_fit = datarun.cones.rf_fits{cell_index};
            
            % skip if it's empty (indicates no fit was done, or the fit collapsed)
            if isempty(the_fit);continue;end
            
            % take center minus surround
            norm_factor = the_fit.center_scale - the_fit.surround_scale;
        
        case {'center sum','center std'} % normalize based on "center" cones, i.e. positive cones within 2 center sigmas
            
            % get the fit
            the_fit = datarun.cones.rf_fits{cell_index};

            % get center point
            rf_ctr = rf_center(datarun,datarun.cell_ids(cell_index),'com');
            
            % skip if no fit or center point
            if isempty(the_fit) || isempty(rf_ctr);continue;end
            
            % get cones within 2 center sigmas
            [x,y] = rf_cone_profile(the_weights,...
                datarun.cones.centers,rf_ctr,'radius',2*the_fit.center_radius);
            
            % get center weights
            center_weights = y(y>0);
            
            
            
            switch params.type
                
                case 'center sum' % normalize to the sum of "center" cones, i.e. positive cones within 2 center sigmas
                    
                    norm_factor = sum(center_weights);
                
                case 'center std' % divide by the std of the center cones
                    
                    norm_factor = std(center_weights);
            
            end
            
        case 'max' % divide by the maximum
            
            % get the max
            norm_factor = max(the_weights);
            
            
        case 'noise std' % divide by the robust std
            
            % get the max
            norm_factor = robust_std(the_weights);
            
            
            
    end
    
    datarun.cones.weights(:,cell_index) = the_weights / norm_factor;
    
end




    
    
    
    
    
    
    