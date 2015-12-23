function datarun=load_map(datarun, varargin) 
% loads map file eg from map_ei_classification_txt and organizes cell_types
%
% datarun{1} master
% datarun{2} slave
% [datarun{..} slave]
%
% default .verbose 1 
%
% greschner

% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addParamValue('field_name','map_path');%
p.addParamValue('verbose', 1);
% parse inputs
p.parse(varargin{:});
params = p.Results;


tmap=dlmread(datarun{2}.names.(params.field_name));


%map cell_ids
    for i=2:size(tmap,2)
        
        datarun{i}.cell_ids_map=zeros(length(datarun{i}.cell_ids),2);
        datarun{i}.cell_ids_map(:,2)=datarun{i}.cell_ids;
        datarun{i}.cell_ids=zeros(length(datarun{i}.cell_ids),1);
        
        for ii=1:size(tmap,1)          
            t=find(datarun{i}.cell_ids_map(:,2)==tmap(ii,i));
            datarun{i}.cell_ids_map(t,1)=tmap(ii,1);
            datarun{i}.cell_ids(t)=tmap(ii,1);            
        end
    end

    
    
%remove unmaped
    for i=1:length(datarun{1}.cell_types)
      
        temp=[];
        for ii=1:length(datarun{1}.cell_types{i}.cell_ids)
            t=find(tmap(:,1)==datarun{1}.cell_types{i}.cell_ids(ii));
            if ~isempty(t) && all(tmap(t,:)) 
                temp=[temp datarun{1}.cell_types{i}.cell_ids(ii)];
            end          
        end
        
        org(i,:)=[length(datarun{1}.cell_types{i}.cell_ids) length(temp)];
        datarun{1}.cell_types{i}.cell_ids=temp;    
    end

  
    
for i=2:size(tmap,2)
    datarun{i}.cell_types=datarun{1}.cell_types;
end  


if params.verbose
    disp('load_map:');
    for i=1:length(datarun{1}.cell_types)
        disp(sprintf('  %3.0f -> %3.0f  %s',org(i,1),org(i,2),datarun{1}.cell_types{i}.name));
    end
end
    
    
 


if 0    
for i=1:length(datarun{1}.cell_types)
    
    for ii=1:size(tmap,2)
    	res{ii}.cell_types{i}.name=datarun{1}.cell_types{i}.name;
        res{ii}.cell_types{i}.cell_ids=[];
    end
    
    for ii=1:length(datarun{1}.cell_types{i}.cell_ids)
        
        t=find(tmap(:,1)==datarun{1}.cell_types{i}.cell_ids(ii));
        
        if all(tmap(t,2:end)) 
            for iii=1:size(tmap,2)
                res{iii}.cell_types{i}.cell_ids=[res{iii}.cell_types{i}.cell_ids tmap(t,iii)];
            end         
        end
 
    end
 
end


if params.verbose
    disp('load_map:');
    for i=1:length(datarun{1}.cell_types)
        disp(sprintf('  %3.0f -> %3.0f  %s',length(datarun{1}.cell_types{i}.cell_ids),length(res{1}.cell_types{i}.cell_ids),datarun{1}.cell_types{i}.name));
    end
end


for i=1:size(tmap,2)
    datarun{i}.cell_types=res{i}.cell_types;
end
end







