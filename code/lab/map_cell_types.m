function datarun=map_cell_types(datarun, varargin) 
% checks which cell_ids from cell_types in master (datarun_1) exist in
% slave (datarun_1) and gives back slave cell_types-structur
%
% .verbose 1
% .map 'All' maps all datarun{..} alternativly array []  
%
% greschner

% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addParamValue('verbose', 1);%
p.addParamValue('map', 'All');%
% parse inputs
p.parse(varargin{:});
params = p.Results;


if isequal(params.map,'All')
    params.map=[1:length(datarun)];
end


for i=1:length(datarun{1}.cell_types)

    cell_types{i}.name=datarun{1}.cell_types{i}.name;

    temp=[];
    for ii=1:length(datarun{1}.cell_types{i}.cell_ids) 
        t=zeros(size(params.map));
        for iii=1:length(params.map)
            t(iii)=ismember(datarun{1}.cell_types{i}.cell_ids(ii),datarun{params.map(iii)}.cell_ids);
        end
        if all(t)            
            temp=[temp datarun{1}.cell_types{i}.cell_ids(ii)];            
        end
    end
    cell_types{i}.cell_ids=temp;

end


if params.verbose
    disp('cell_type mapping:');
    for i=1:length(datarun{1}.cell_types)
        disp(sprintf('  %3.0f -> %3.0f  %s',length(datarun{1}.cell_types{i}.cell_ids),length(cell_types{i}.cell_ids),datarun{1}.cell_types{i}.name));
    end
end


for i=1:length(params.map)
    datarun{params.map(i)}.cell_types=cell_types;
end



if 0
    if length(datarun)==2 | isequal(params.map,[1 2])
        for i=1:length(datarun{1}.cell_types)

            cell_types{i}.name=datarun{1}.cell_types{i}.name;

            temp=[];
            for ii=1:length(datarun{1}.cell_types{i}.cell_ids)        
                if ismember(datarun{1}.cell_types{i}.cell_ids(ii),datarun{2}.cell_ids)            
                    temp=[temp datarun{1}.cell_types{i}.cell_ids(ii)];            
                end
            end
            cell_types{i}.cell_ids=temp;

        end
    end     
end














