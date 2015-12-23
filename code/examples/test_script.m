
% This script contains some examples of the proposed functions to show how they work


% load and look at some STAs
if 1
    % initialize datarun
    datarun = load_data('2005-08-08-0','rf-1-auto-jg-0');
    
    % load params file
    datarun = load_params(datarun,'verbose',1);
    
    % load sta file
    datarun = load_sta(datarun,'load_sta',[],'verbose',1,'keep_java_sta',1);
    
    % load the various STA summaries, don't save STAs
    datarun = get_sta_summaries(datarun,{1},'verbose',1,'keep_stas',0,'keep_rfs',1);
    
    % set polarities
    datarun = set_polarities(datarun);
    
    % get a quick look at ON parasol cells
    plot_rf_portraits(datarun,{1},'plot_radius',5);
    
    % load vision fits
    datarun = get_sta_fits_from_vision(datarun,'all');
    
    % plot ON parasol fits 
    figure;plot_rf_summaries(datarun,{1},'plot_fits',1,'label',1)
end


% compute and plot SNLs
if 0
    % load datarun
    datarun = load_data('2005-08-08-0/data010-4dot5/data010-4dot5');
    
    % load STAs, classification, and spikes
    datarun = load_params(datarun,'verbose',1);
    datarun = load_sta(datarun,'verbose',1,'load_sta',{1});
    datarun = load_neurons(datarun);
    datarun = get_sta_summaries(datarun,{1},'verbose',1,'fig_or_axes',1);
    
    % create movie
    datarun = load_java_movie(datarun,'/snle/acquisition/movie-xml/RGB-10-8-0.48-33333.xml');
    
    % compute SNLs
    datarun = get_snls(datarun,{1},'end_stim',1000);
    
    % plot SNLs
    plot_snls(datarun,{1})
end



% load series information
%   file names
%   stimulus
%   piece info

% load neurons file
%   spikes
%   other stuff

% load params file
%   vision classification
%   STA fits
%   other stuff

% load stas file

% load eis file





% BORING MECHANICAL STUFF

% test case handling of load_data
if 0
    % vision path
    datarun = load_data('/vison/path');datarun.names
    datarun = load_data('/vison/path',struct('verbose',true));datarun.names
    % expt & condition
    datarun = load_data('expt','cond');datarun.names
    datarun = load_data('expt','cond',struct('verbose',true));datarun.names
    % struct
    datarun = load_data(struct('names',struct('experiment','expt','condition','cond')));datarun.names
    datarun = load_data(struct('names',struct('experiment','expt','condition','cond')),struct('verbose',true));datarun.names
    % multiple specifications
    datarun = load_data({{'expt1','cond1'},{'expt2','cond2'}});datarun{1}.names,datarun{2}.names
    datarun = load_data({{'expt1','cond1'},{'expt2','cond2'}},struct('verbose',true));datarun{1}.names,datarun{2}.names
    
end

% test argument handling of default_params
if 0
    defaults = struct('a',1,'b',2,'c',3);
    
    % empty params struct
    default_params(defaults,[])
    
    % specify a value
    default_params(defaults,struct('b',2.2))
    
    % specify a value incorrectly
    default_params(defaults,struct('bb',2.2))
    
end

% get_cell_indices
if 0
    ct{1} = struct('name','ON parasol','cell_ids',datarun.cell_ids([1 3 5]));
    ct{2} = struct('name','OFF parasol','cell_ids',datarun.cell_ids([2 4 6]));
    datarun.cell_types = ct;
    
    get_cell_indices( datarun, 5851 )
    get_cell_indices( datarun, [5851 6588] )
    get_cell_indices( datarun, {1} )
    get_cell_indices( datarun, 'ON parasol' )
    a=get_cell_indices( datarun, 'ON parasol' )
    [a,b]=get_cell_indices( datarun, 'ON parasol' )
    [a,b,c]=get_cell_indices( datarun, 'ON parasol' )
    [a,b,c] = get_cell_indices( datarun, {1,2} )
    [a,b,c] = get_cell_indices( datarun, {'ON parasol',2} )

end
    

