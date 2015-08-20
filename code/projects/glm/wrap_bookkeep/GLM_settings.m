function GLMType = GLM_settings(base_type, changes_cell)
%%% PURPOSE %%%
% Save space in other codes for dictating GLM Parameters and names
% Code which dictates the GLMType name which gets interpreted downstream
% All names for CLMTypes should come from here.


%%% INPUTS  %%%
% Base_Type
%  for now just use default to this! 
% Changes_cell
%

%%% OUTPUTS %%%
% GLMType
% Structure which uniquely defines GLM Type (up to stim_type)
% .map_type: mapping from master (classification) to slave
%            refers to name associated with actual raw data structures in vision                    
% .cone_model : type of input filtering used before the GLM
% .cone_sname : shorter version used to identify input model 
%               sname stands for "short name"
% .nullpoint  : where we set our input to be zero .. 'mean'
%               so far is only set to 'mean'  could change later!
% .stimfilter_mode: dictates GLM spatial temporal stimulus type
%                  ie. rk-1 , fixedSP_rk1_linear,  etc.
% .CONVEX : a direct consequence of .stimfilter_mode
%           true or false statement,  rk1 false, 
% .TonicDrive
% .StimFilter
% .PostSpikeFilter
% .nullpoint

%%% CALLS %%%
% none


% AKHEITMAN 2014-12-13

if strcmp(base_type, 'default')
    GLMType.stimfilter_mode = 'fixedSP_rk1_linear'; GLMType.CONVEX = true;
    GLMType.TonicDrive = true;
    GLMType.StimFilter = true;
    GLMType.PostSpikeFilter = true;
    GLMType.CouplingFilters = false;
    GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';
    GLMType.nullpoint  = 'mean'; 
    GLMType.map_type   = 'mapPRJ'; 
    GLMType.debug      = false;
    GLMType.contrast   = false;
    GLMType.Subunits   = false;
end
%%%%% Cone Names %%%%%%%

%'timekernel';




if exist('changes_cell','var') && length(changes_cell)>=1
    for i_change = 1:length(changes_cell)
        change = changes_cell{i_change};
        
        if strcmp(change.type, 'cone_model')
            if strcmp(change.name, 'rieke_linear')
                GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
            end
            
            if strcmp(change.name, 'rieke_fullcone')
                GLMType.conemodel = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
            end
            %GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
            %GLMType.cone_model = 'linear_timekernel_shift0';  GLMType.cone_sname ='timekernel';
            %GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
            %GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';
            %GLMType.cone_model = 'linear_timekernel_shift1';  GLMType.cone_sname =
            %GLMType.cone_model = 'TimeConvolve_DimFlash_shift2';
            %GLMType.cone_sname = 'timesmoothed';
        end
        
        if strcmp(change.type, 'filter_mode')
            if strcmp(change.name, 'nostim')
                GLMType.stimfilter_mode = 'nostim'; 
            end
            
            if strcmp(change.name, 'rk2')
                GLMType.stimfilter_mode = 'rk2'; 
                GLMType.CONVEX = false;
            end
            
            if strcmp(change.name, 'rk1')
                GLMType.stimfilter_mode = 'rk1'; 
                GLMType.CONVEX = false;
            end
            if strcmp(change.name,'fixedSP-ConductanceBased')
                GLMType.stimfilter_mode = 'fixedSP-ConductanceBased';
                GLMType.CONVEX = true;
                GLMType.postfilter_nonlinearity      =  true;
                GLMType.postfilter_nonlinearity_type ='ConductanceBased_HardRect';
            end
            if strcmp(change.name, 'rk2-ConductanceBased')
                GLMType.stimfilter_mode = 'rk2-ConductanceBased'; 
                GLMType.CONVEX = false;
                GLMType.postfilter_nonlinearity      =  true;
                GLMType.postfilter_nonlinearity_type ='ConductanceBased_HardRect';
            end
            %GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; 
            %GLMType.fixedSPlength = 13;  
            %GLMType.fixedSP_nullpoint = 'mean'; 
            %GLMType.stimfilter_mode = 'rk2-ConductanceBased';
            %GLMType.stimfilter_mode = 'rk2';
            %GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
        end

        
         if strcmp(change.type, 'postfilter_nonlinearity')
            GLMType.postfilter_nonlinearity      =  true;
            if strcmp(change.name, 'piece_linear_aboutzero')
                GLMType.postfilter_nonlinearity_type =  'piece_linear_aboutzero';
                GLMType.DoubleOpt = true;
                GLMType.DoubleOpt_Manual = true;
            end
         end
        

        if strcmp(change.type, 'debug')
            if strcmp(change.name, 'true')
                GLMType.debug = true;
            elseif strcmp(change.name,'false')
                GLMType.debug = false;
            else
                error('wrong name for debug setting')
            end
        end
        
        %
        if strcmp(change.type, 'input_pt_nonlinearity')
            if strcmp(change.name, 'piece_linear_aboutmean')
                 GLMType.input_pt_nonlinearity        = true;
                 GLMType.input_pt_nonlinearity_type   = 'piece_linear_aboutmean';
                 GLMType.DoubleOpt = true;
                 GLMType.DoubleOpt_Manual = true;
            end
            if strcmp(change.name, 'piecelinear_fourpiece_eightlevels')
                 GLMType.input_pt_nonlinearity        = true;
                 GLMType.input_pt_nonlinearity_type   = 'piecelinear_fourpiece_eightlevels';
                 GLMType.DoubleOpt = true;
                 GLMType.DoubleOpt_Manual = true;
            end
            %  AKH 2015-07-14  added power raise option
            if strcmp(change.name, 'log_powerraise')
                 GLMType.input_pt_nonlinearity        = true;
                 GLMType.input_pt_nonlinearity_type   = 'log_powerraise';
                 GLMType.InputNL_IteratedOpt = true;
            end
        end
        if strcmp(change.type, 'PostSpikeFilter')
            if strcmp(change.name, 'OFF')
                GLMType.PostSpikeFilter = false;
            end
        end
        
        if strcmp(change.type, 'CouplingFilters')
            if strcmp(change.name, 'ON')
                GLMType.CouplingFilters = true;
            end
        end
        
        if strcmp(change.type, 'specialchange')
            GLMType.specialchange = true;
            GLMType.specialchange_name = change.name;
        end
        
        if strcmp(change.type, 'Subunits') && strcmp(change.name, 'ON')
            GLMType.Subunits = true;
        end
        
        if strcmp(change.type, 'Contrast') && strcmp(change.name, 'ON')
            GLMType.contrast = true;
        end

        %{
        %GLMType.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
        %GLMType.input_pt_nonlinearity_type = 'piece_linear_shiftmean';
        %GLMType.input_pt_nonlinearity_type = 'polynomial_order5_part4';
        %GLMType.input_pt_nonlinearity_type = 'piecelinear_fourpiece_eightlevels';
        %GLMType.input_pt_nonlinearity_type =  'oddfunc_powerraise_aboutmean';
        %GLMType.input_pt_nonlinearity_type =  'log';
        %GLMType.input_pt_nonlinearity_type =  'exp';
        %GLMType.input_pt_nonlinearity_type = 'polynomial_androot_order2_search2';  % order plus minus 2
        %GLMType.input_pt_nonlinearity_type = 'polynomial_androot_order2_search3';

        %GLMType.postfilter_nonlinearity      =  true;
        %GLMType.postfilter_nonlinearity_type = 'ConductanceBased_HardRect';
        %GLMType.postfilter_nonlinearity_type =  'oddfunc_powerraise_aboutmean';
        %GLMType.postfilter_nonlinearity_type =  'piece_linear_aboutmean';
        %GLMType.postfilter_nonlinearity_type = 'raisepower_meanafter';
        %Type = 'Stim_Nonlinearity'; modification = 'raisepower_meanafter'
        %GLMType.CONVEX = false; % with relation to the filters .. are parameters used linearly in the GLM. 
        %GLMType.DoubleOpt = true;
        %GLMType.DoubleOpt_Manual = true;
        %GLMType.stimfilter_mode = 'rk1';
        %GLMType.specialchange = false;
        %GLMType.specialchange_name = 'Conductance_Based';
        %}

        
        %GLMType.fixed_spatialfilter = true;
    end
end
       


end
