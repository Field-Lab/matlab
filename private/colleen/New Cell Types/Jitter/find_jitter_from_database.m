% find jitter pieces

[~,~,RAW_data]=xlsread('/Volumes/Lab/Users/crhoades/Database spreadsheet.xlsx', 'Dataruns');
jitter =[];
for i = 1:size(RAW_data,1)
    if strfind(RAW_data{i,31}, 't')
        jitter = [jitter; i];
    end
end


for i = 1:length(jitter)
    test_index = jitter(i);
    date{i,2} = RAW_data{test_index,19};
    date{i,3} = RAW_data{test_index,13};

    while true
        
        if ~isnan(RAW_data{test_index,1})
            date{i,1} = RAW_data{test_index,1};

            break;
        else
            test_index = test_index-1;
        end
        
    end
end

    