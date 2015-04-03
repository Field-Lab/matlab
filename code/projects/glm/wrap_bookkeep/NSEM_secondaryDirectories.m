% sample
%{
directory_type = 'GLM_fitandraster';
inputs.exp_nm = '2012-08-09-3';
inputs.stim_type = 'WN';
inputs.map_type = 'mapPRJ';
condition = 'alligator';
directory = NSEM_secondaryDirectories(directory_type, inputs, condition)
%}

function directory = NSEM_secondaryDirectories(directory_type, inputs, condition,BaseDirectories)
% 2014-04-05 AKHeitman

% Condtion is an optional argument meant for transitory perioud between
% Salk/ Stanford

% Condition can be either netappsnle  or alligator

%%%%%inputs
%.GLMoutput  -GLM output directory, should be dictated by NSEM_BaseDirectories
%            -default BD.GLM_output;
%.exp_nm  - 
%.fit_type   - either 'NSEM' or 'BW'  ('or WN')



if exist('BaseDirectories', 'var')
    BD = BaseDirectories;
else
    BD = NSEM_BaseDirectories;
end

GLM_Output_Raw       = BD.GLM_output_raw;
GLM_Output_Analysis  = BD.GLM_output_analysis;

BlockedSpikes        = BD.BlockedSpikes;

% temporary directory discrepancy %
if exist('condition','var')
    if strcmp(condition, 'netappsnle')
        GLM_Output_Raw = BD.temp.GLM_output_netappsnle;
    end
    if strcmp(condition, 'alligator')
        GLM_Output_Raw = BD.temp.GLM_output_alligator;
    end
end

% how to account for old conditions where WN was called BW
if isfield(inputs, 'stim_type') 
    if exist('condition','var') && strcmp(condition, 'netappsnle')
        if strcmp(inputs.stim_type, 'WN')
            inputs.stim_type = 'BW';
        end
    end
end



switch directory_type
	case 'loaddir_GLMfit'
        % location for Loading GLMFits for ANalysis
        % GLMFits which have been deemed accurate (rather than raw fits)
        % Run out of the GLM_OUtput_Analysis
        % Builds off of BD.GLM_output
        % Dictated by
        %    -piece number
        %    -the fitting stimulus
        %    -how the master maps to the slave
        %    -specific GLM fit params
        
        exp_nm          = inputs.exp_nm;
        stim_type       = inputs.stim_type;
        map_type        = inputs.map_type;
        fitname         = inputs.fitname;
        directory = sprintf('%s/%s/%s_%s/%s', GLM_Output_Analysis, fitname, stim_type, map_type,exp_nm);
        %{
        if exist('condition', 'var')
            GLM_fitandraster       = NSEM_secondaryDirectories('GLM_fitandraster', inputs, condition);
        else
            GLM_fitandraster       = NSEM_secondaryDirectories('GLM_fitandraster', inputs);
        end
        
        directory = sprintf('%s/%s', GLM_fitandraster , paramname);
        %}
    
    case 'savedir_GLMfit'
        % location for saving individual fits. 
        % Builds off of BD.GLM_output
        % Dictated by
        %    -piece number
        %    -the fitting stimulus
        %    -how the master maps to the slave
        %    -specific GLM fit params
        
        exp_nm          = inputs.exp_nm;
        stim_type       = inputs.stim_type;
        map_type        = inputs.map_type;
        fitname         = inputs.fitname;
        directory = sprintf('%s/%s/%s_%s/%s', GLM_Output_Raw, fitname, stim_type, map_type,exp_nm);
        %{
        if exist('condition', 'var')
            GLM_fitandraster       = NSEM_secondaryDirectories('GLM_fitandraster', inputs, condition);
        else
            GLM_fitandraster       = NSEM_secondaryDirectories('GLM_fitandraster', inputs);
        end
        
        directory = sprintf('%s/%s', GLM_fitandraster , paramname);
        %}
    case 'organizedspikes_dir'
        % location for saving individual fits. 
        % Builds off of BD.GLM_output
        % Dictated by
        %    -piece number
        %    -the fitting stimulus
        %    -how the master maps to the slave
        %    -specific GLM fit params

        exp_nm          = inputs.exp_nm;
        stim_type       = inputs.stim_type;
        map_type        = inputs.map_type;
        directory = sprintf('%s/%s/%s_%s',BlockedSpikes, exp_nm, stim_type, map_type);
    
    case 'WN_STA'
        exp_nm          = inputs.exp_nm;
        map_type        = inputs.map_type;
        directory = sprintf('%s/%s/%s_%s/STA',BlockedSpikes, exp_nm, 'WN', map_type);
        
 %}
        
end
        
