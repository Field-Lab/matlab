% loop through text files on disk, update things (e.g. change cone classification),
% and write files with the new info


% add java path
javaaddpath /snle/lab/Applications/Vision.app/Contents/Resources/Java/Vision.jar
clear java


% IDENTIFY SAVED DATARUNS 

% get list of files
dir_contents = dir(single_cone_path);
files_to_load = cell(0);
% look though each file
for ff=1:length(dir_contents)
    % if it's a saved datarun (i.e. starts with "save-")
    if strfind(dir_contents(ff).name,'save-') == 1
        % check for presence of text files
        short_name = strrep(strrep(dir_contents(ff).name,'save-',''),'.mat','');
        if exist([single_cone_path short_name '/cones.txt'],'file')
            % if they exist, load this mat file
            files_to_load{length(files_to_load)+1} = dir_contents(ff).name;
        end
    end
end
clear dir_contents ff short_name
% note the quantity
fprintf('\nWill load %d dataruns\n\n',length(files_to_load))



% LOOP THROUGH THEM

for dd = 1:length(files_to_load)

    fprintf('\n ****************   Now loading %s  **************** \n',files_to_load{dd})

    try

        % load datarun
        load([single_cone_path files_to_load{dd}])

        % ensure proper fields exist (by filling them in with BS values!)
        if ~isfield(datarun.piece,'optical_path_direction');datarun.piece.optical_path_direction = '?';end
        if ~isfield(datarun.piece,'rig');datarun.piece.rig = '?';end

        % load cone info
        [datarun,extra_info] = import_single_cone_data(datarun);

        % get short name
        short_name = strrep(strrep(files_to_load{dd},'save-',''),'.mat','');

        % make changes
        switch 1
            case 1  % reclassify cones, make cone ROI
                datarun.stas.rfs = [];

                % parameters
                num_bins = 20;
                bin_width = 1;
                save_file_path = [single_cone_path short_name '/'];

                % classify cones
                [cone_types,likelihood_ratios,classify_extras] = ...
                    classify_cones(datarun.cones.rgb, cone_rgb_expected(datarun),'algorithm','k means, EM');
                
                % enter classification into datarun
                datarun.cones.types = cone_types;
                datarun.cones.likelihoods = likelihood_ratios;
                datarun.cones.types_em = classify_extras.types_em;
                datarun.cones.types_kmeans = classify_extras.types_kmeans;

                % get local density
                [drp,bin_centers,extras] = density_recovery_profile(datarun.cones.centers,num_bins,bin_width);

                % identify ROI
                cone_roi = identify_cone_mosaic_roi(datarun.cones.centers,extras.density, extras.eff_rad);

                % plot cone mosaics and classification figures
                plot_cone_summary(datarun,'fig_drp',1,'fig_class',2,'cone_roi',cone_roi);

                % print figures
                print(1,[save_file_path 'mosaic'],'-dpdf')
                print(2,[save_file_path 'classification'],'-dpdf')
        end

        % if the ROIs were not assigned above, assign them to whatever was read from disk
        if ~exist('cone_roi','var');cone_roi = extra_info.cone_roi; end
        if ~exist('rgc_roi','var');rgc_roi = extra_info.rgc_roi; end

        % write to disk
        export_single_cone_data(datarun,extra_info.cell_ids,[],save_file_path,'rgc_roi',rgc_roi,'cone_roi',cone_roi)

    catch

        % get the error
        temp=lasterror;

        % display the error and some eye-catching text
        fprintf('\n\n\n******* ERROR in %s *******\n\n',files_to_load{dd})
        disp(temp.message)
        disp(temp.identifier)
        for ss =1:length(temp.stack)
            disp(temp.stack(ss))
        end
        fprintf('\n\n\n')

    end

    % clean up
    clear datarun short_name cone_roi rgc_roi extras extra_info save_file_path drp

end



        
% remove java path
javarmpath /snle/lab/Applications/Vision.app/Contents/Resources/Java/Vision.jar
clear java

