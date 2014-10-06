function save_out_jitter_stas(datarun, cell_spec, save_dir, params)
% SAVE_OUT_JITTER_STAS     For many cells, compute and save the STA for each stimulus offset
%
% usage:  save_out_jitter_stas(datarun, cell_spec, save_dir, params)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells
%              params - struct of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional fields in params, their default values, and what they specify:
%
% summary           true	also compute and save the jittered STA summary frame
% frames            -3:0    which frames to save out
% spike_fraction   	0.8     what fraction of spikes to use
% iterations      	1    number of times to compute the STA for each cell
%
%
% example usage:
%
%   save_dir = ['~gauthier/Desktop/2008-one/code/jitter_' datarun.names.datarun_name];
%   save_out_jitter_stas(datarun, 1621, save_dir, struct('spike_fraction',0.8,'iterations',1000))
%
%
% 2008-09  gauthier
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.summary = true;
defaults.frames = -3:0;
defaults.spike_fraction = 1;
defaults.iterations = 1;

% combine user and default parameters
params = default_params( defaults, params);




% BODY OF THE FUNCTION


% define cells to loop through
cell_nums = get_cell_indices(datarun,cell_spec);

% do more iterations if desired
cell_nums = repmat(cell_nums,1,params.iterations);

% get offset_index
offset_index = datarun.stas.jitter.offset_index;


% generate white noise movie
fprintf('\nCalculating white noise movie...')
start_loading = clock; % note when it started
import('edu.salk.snl.cellfinder.statistics.*');
bufMov = ReceptiveField.calculateMovie(datarun.stimulus.mdf_file,int32(datarun.triggers*datarun.sampling_rate));

% note how long it took
fprintf(' done (%0.1f minutes)\n\n',etime(clock,start_loading)/60);



% loop through cells
for cc = 1:length(cell_nums)

    %try

    % get cell id
    cell_id = datarun.cell_ids(cell_nums(cc));

    % reset buffered movie
    bufMov.reset

    
    
    % compute STA for each jitter
    
    % either using all spikes
    if params.spike_fraction == 1 
        jitter_stas = compute_jitter_sta(datarun,offset_index,cell_id,bufMov,params.frames);
        
        
    else % or a subset
        
        % set random seed
        rand_seed = round(mod(sum(clock*100),10000));
        
        % display it
        fprintf('\nrandom seed is %d\n',rand_seed)
        
        % compute the jittered STAs
        jitter_stas = compute_jitter_sta(datarun,offset_index,cell_id,bufMov,params.frames,...
            struct('spike_fraction',params.spike_fraction,'rand_seed',round(mod(sum(clock*100),10000))));
    end

    

    % for each jittered STA, compute summary frame, if desired
    
    fprintf('#####\n')
    keyboard;
    if params.summary

        % get timecourse from coarse sta
        sta_coarse = datarun.stas.stas{get_cell_indices(datarun,cell_id)};
        tc =  time_course_from_sta(sta_coarse,struct('polarity',datarun.stas.polarities(cell_nums(cc))));

        % use only the frames that are kept from the jittered STAs
        tc = tc(size(tc,1)+params.frames,:);

        % regress against this timecourse for each sta
        for oo = 1:size(offset_index,2)
            jitter_summaries{oo} = rf_from_sta(jitter_stas{oo},'regress',struct('r_timecourse',tc));

            % keep only green (FOR NOW!)
            %jitter_summaries{oo} = jitter_summaries{oo}(:,:,2);
        end
    end


    
    % save to disk

    % set name of directory
    if params.spike_fraction == 1 
        stas_dir = save_dir;
    else
        stas_dir = sprintf('%s/%d-subset-%s/',save_dir,cell_id,strrep(sprintf('%0.3f',params.spike_fraction),'.','_'));
    end
    
    % ensure the directory exists
    if ~exist(stas_dir); mkdir(stas_dir); end
    
    % set the name where to save
    if params.spike_fraction == 1 
        save_name = sprintf('%s/%d.mat',stas_dir,cell_id);
    else
        save_name = sprintf('%s/%d-%d.mat',stas_dir,cell_id,rand_seed);
    end
    
    % save stas and summary, if desired
    if params.summary
        save(save_name,'jitter_stas','jitter_summaries');
    else
        save(save_name,'jitter_stas');
    end
        
    %catch
    %    disp(lasterr)
    %end
end


