
    
    %%
    data1 = '2015-04-14-2';
    run1 = 'data008';
     datarun.names.rrs_params_path=['/Volumes/Acquisition/Analysis/', data1, '/', run1,'/',run1, '.params'];
 datarun.names.rrs_neurons_path=['/Volumes/Acquisition/Analysis/', data1, '/', run1, '/',run1,'.neurons'];
  datarun.names.rrs_ei_path=['/Volumes/Acquisition/Analysis/', data1, '/', run1, '/',run1,'.ei'];

 opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_ei', true);

    datarun1= load_data(datarun,opt);
    
        data2 = '2015-04-14-2';
    run2 = 'data008';
     datarun.names.rrs_params_path=['/Volumes/Analysis/', data2, '/', run2, '/',run2,'.params'];
 datarun.names.rrs_neurons_path=['/Volumes/Analysis/', data2, '/', run2,'/', run2,'.neurons'];
  datarun.names.rrs_ei_path=['/Volumes/Analysis/', data2, '/', run2,'/', run2,'.ei'];
 opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_ei', true);

    datarun2= load_data(datarun,opt);
    
    [cell_list_map, failed_cells, output_matrix] = map_ei(datarun2, datarun1, 'corr_threshold', 0.8, 'electrode_threshold', 4);
    
    master_class = importdata(['/Volumes/Analysis/', data1,'/', run1, '/', run1, '.classification.txt']);
    for i = 1:size(master_class,1)
        row = master_class{i};
        s = strfind(row, ' ');
        cell_id(i,1) = str2num(row(1:(s(1)-1)));
        category{i,1} = row(s(end)+1:end);
    end
    i = 0
    filename = [data1, ' ', run2, ' EIs from ', run1,'.txt']; 
    fid = fopen(filename, 'a');
    for i=1:size(output_matrix,1)
        if output_matrix(i,1) ~= 0
            ind = cell_id == output_matrix(i,2);
            cat = category{ind};
            toAddtoFile = [num2str(output_matrix(i,1)), '  ', cat]; 
            fprintf(fid, [toAddtoFile, '\n']);
%             ind = cell_id == output_matrix(i,3);
%             cat = 'All/duplicates/'
            i = i+1;
        end
    end
    
        fclose(fid)
        
            
               
