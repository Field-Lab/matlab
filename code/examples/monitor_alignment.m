% align the monitor with the array
%
%
% detailed documentation of photographic mapping in
%
%   matlab-standard/documentation/monitor to array alignment.intaglio
%
%
%
% 2010-04  gauthier
%
%



%%
%%%%%%%% PHOTOGRAPHIC MAPPING %%%%%%%%%


% choose dataset

%datarun = load_data('2007-09-18-4/data002-nwpca/data002');      %      512 recording
%datarun = load_data('2009-04-13-1/data013/data013/data013');    %      519, 30 um, 10x from below
datarun = load_data('2010-03-05-2/data000/data000');            %      519, 30 �m, 6.5x from below

% load datarun
datarun = load_index(datarun);
datarun = load_sta(datarun,'load_sta',[]);
datarun = load_params(datarun);
datarun = load_ei(datarun,[]);
datarun = get_sta_fits_from_vision(datarun);

% load alignment information from disk
datarun = load_monitor_alignment(datarun);


% OPTIONAL: adjust the alignment points
if 0
    datarun = compute_monitor_to_array_transformation(datarun);
    % save the result
    %save_monitor_alignment(datarun);
end

% verify the alignment points look good
verify_array_monitor_alignment(datarun)


% plot many RFs with the array outline
figure;plot_rf_summaries(datarun,'all','label',true,'array',true,'array_label',true)

% plot an RF with the array
plot_rf(datarun,datarun.cell_ids(1),'array',true)
% add the EI
hold on;plot_ei(datarun,datarun.cell_ids(1),'coordinates','sta','pretty_axes',0,'alpha',0)

% plot an STA
plot_sta(datarun,datarun.cell_ids(1),'array',true)

% plot many RFs with their EIs
plot_rfs_eis(datarun,{1})





%%
%%%%%%%% CLICKED POINT MAPPING %%%%%%%%%


% choose dataset

datarun = load_data('2008-08-27-2/data001/data001');   %               61, Rig C


% load datarun
datarun = load_index(datarun);
datarun = load_sta(datarun,'load_sta',[]);
datarun = load_params(datarun);
datarun = load_ei(datarun,[]);
datarun = get_sta_fits_from_vision(datarun);


% compute the transformation using the clicked points in datarun.piece.array
datarun = compute_monitor_to_array_transformation(datarun);


% verify the interpolation of clicked points look good
verify_array_monitor_alignment(datarun)



% plot many RFs with the array outline
figure;plot_rf_summaries(datarun,{3},'label',true,'array',true,'array_label',true)

% plot an RF with the array
plot_rf(datarun,datarun.cell_ids(1),'array',true)
% add the EI
hold on;plot_ei(datarun,datarun.cell_ids(1),'coordinates','sta','pretty_axes',0,'alpha',0)

% plot an STA
plot_sta(datarun,datarun.cell_ids(1),'array',true)

% plot many RFs with their EIs
plot_rfs_eis(datarun,{1})



