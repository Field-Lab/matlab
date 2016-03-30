load area_new.mat

% [area] = compare_RF_sizes(1);
[match] = get_cell_type_index(34);


[~, txt,array_size] = xlsread('/Volumes/Lab/Users/crhoades/Large Cell Data ARVO.xlsx');
array_size = cell2mat(array_size(2:end,4));   

pieces = [5 6 17 20 24 27 31 34]

% Compare directly with OFF smooth
% pieces = [9];

for j= 1:length(pieces)
    % piece = txt(j,1:3);
    run_opts.date=strtrim(txt{pieces(j),1}); % one slash at the end
    temp1 = strtrim(txt{pieces(j),2});
    temp2 =  strtrim(txt{pieces(j),3});
    run_opts.concatname=temp1; % Name (or modified name) of run, no slashes\
    run_opts.dataname = temp2;
    run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',run_opts.dataname, '/',run_opts.dataname];
    
    run_opts.save_location_root = '/Volumes/Lab/Users/crhoades/Cell Properties/';
    run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, '/', run_opts.dataname];
    
    output_large = load([run_opts.filepath, '/' 'ON large 1','/output.mat']);
    
    output_parasol= load([run_opts.filepath, '/' 'ON parasol','/output.mat']);
    output_midget= load([run_opts.filepath, '/' 'ON midget','/output.mat']);
    
    parameters_large{pieces(j)}  =output_large.output.parameters;
    parameters_parasol{pieces(j)}  =output_parasol.output.parameters;
    parameters_midget{pieces(j)}  =output_midget.output.parameters;
    
   
    
    parasol_rf(1,j) = mean(area{pieces(j)}{match(pieces(j),1)});
    
    large_rf(1,j) = mean(area{pieces(j)}{match(pieces(j),5)});
    midget_rf(1,j) = mean(area{pieces(j)}{match(pieces(j),3)});
    



[X,Y] = pol2cart(all_pieces{i,3}+pi/2, all_pieces{i,2}); % referenc epoint of all_pieces{i,3} is 12:00, change to be 3:00 by adding pi/2
        test = [test;all_pieces{i,3}];
        if all_pieces{i,3} > (pi)
            E = sqrt((X/0.61)^2+Y^2); %nasal

        else

            E = sqrt(X^2+Y^2); %temporal
            
        end
        %     all_pieces{i,4} = 0.1+4.21*E+0.038*E^2; %visual angle
        all_pieces{i,4} = E; % converted ecc
end
