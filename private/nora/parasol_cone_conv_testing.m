thresh=0.2;

single_cone_datasets={...
'2011-05-11-6',...
'2011-06-24-6',...
'2011-10-25-5',...
'2011-12-13-2',...
'2012-04-13-1',...
'2012-08-21-2',...
'2012-09-06-0',...
'2012-09-13-2',...
'2012-09-21-2',...
'2012-09-24-5',...
'2013-08-19-4',...
'2013-08-19-5 ',...
'2013-10-10-1',...
'2014-04-10-2',...
'2014-06-04-7'};

single_cone_eccentricity=[];

for i_dataset=1:1 %length(single_cone_datasets)
    
    % Find cone data
    cone_data=dir(['/Volumes/Analysis/' single_cone_datasets{i_dataset} '/cone_data']);
    for i_conerun=3:3 %length(cone_data)
        found=false; i=1;
        while ~found
            if strcmp(cone_data(i_conerun).name(i:i+3),'data')
                piece=cone_data(i_conerun).name(i:i+6);
                found=true;
            else
                i=i+1;
            end
        end
        datarun=load_data([single_cone_datasets{i_dataset} '/' piece]);
        datarun=load_neurons(datarun);
        datarun=load_params(datarun);
        datarun=load_cones_njb(datarun, cone_data(i_conerun).name);
        
        % find on and off parasols
        if 0
        for i_type=1:length(datarun.cell_types)
            %$ datarun.cell_types{i_type}.name
            if any(strncmpi(datarun.cell_types{i_type}.name, {'onpar', 'on par'},5))
                on_par=datarun.cell_types{i_type}.cell_ids;
                disp(datarun.cell_types{i_type}.name)
            elseif any(strncmpi(datarun.cell_types{i_type}.name, {'offpar', 'off par'},5))
                off_par=datarun.cell_types{i_type}.cell_ids;
                disp(datarun.cell_types{i_type}.name)
            end
        end
        end
        
        [indices,~,~ ]=get_cell_indices(datarun, {'ON parasol', 'OFF parasol'});
        
        ON_cones=zeros(length(indices{1}),1);
        OFF_cones=zeros(length(indices{2}),1);
        
        % find # of cones for each
        for i_cell=1:length(indices{1})
            cone_weights=datarun.cones.weights(:,indices{1}(i_cell));
            ON_cones(i_cell)=sum(cone_weights>max(cone_weights)*thresh);
            clear cone_weights
        end
        for i_cell=1:length(indices{2})
            cone_weights=datarun.cones.weights(:,indices{2}(i_cell));
            OFF_cones(i_cell)=sum(cone_weights>max(cone_weights)*thresh);
            clear cone_weights
        end
        
        % save in some data structure
        cone_data{i_dataset}.ON.mean=mean(ON_cones);
        cone_data{i_dataset}.ON.stderror=std(ON_cones)/sqrt(length(ON_cones));
        cone_data{i_dataset}.OFF.mean=mean(OFF_cones);
        cone_data{i_dataset}.OFF.stderror=std(OFF_cones)/sqrt(length(OFF_cones));
        
    end
end