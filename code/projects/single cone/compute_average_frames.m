% loop through dataruns on disk, compute average stimulus frame and save it


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
        
        % get short name
        short_name = strrep(strrep(files_to_load{dd},'save-',''),'.mat','');
        
        % get path to files
        rrs_prefix = strrep(short_name,'_','/');
        
        
        if ~strcmp(rrs_prefix,...
                '2007-03-27-2/data014/data014/data014')
            continue
        end
        
        % load datarun based on this
        datarun = load_data(rrs_prefix);
        
        
        % set save path
        save_file_path = [single_cone_path short_name '/'];
        
        switch 2
            case 1 % call vision to compute average frame
                % compute average movie frame
                frame = compute_average_white_noise_frame(datarun);
                
            case 2 % average all STAs together!
                
                % load stas object
                datarun = load_sta(datarun,'verbose',1,'load_sta',[]);
                
                % go through each cell
                cell_indices = get_cell_indices(datarun,'all');
                for cc = 1:length(cell_indices)
                    
                    fprintf('.')
                    
                    % get sta
                    sta=get_sta(datarun,datarun.cell_ids(cell_indices(cc)));
                    
                    % skip if empty
                    if isempty(sta);continue;end;

                    %normlize noise sigma
                    noise_sigmas = robust_std(reshape(permute(sta,[1 2 4 3]),[],size(sta,3)));
                    for ll =1:size(sta,3)
                        sta(:,:,ll,:) = sta(:,:,ll,:)/noise_sigmas(ll);
                    end

                    % get first frame
                    sta = sta(:,:,:,1);

                    if ~exist('sta_sum','var')
                        sta_sum = sta;
                    else
                        sta_sum = sta_sum + sta;
                    end
                end
                
                fprintf('\n')

                sta_avg = sta_sum/length(cell_indices);

        end

        % save to disk
        save([save_file_path 'average_sta_frame'],'sta_avg')
      
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

