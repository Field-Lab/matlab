thresh=0.4;
%

single_cone_datasets={...
'2012-09-21-2/data007',...
'2011-05-11-6/data002',...
'2012-09-24-5/data001',...
'2011-06-24-6/data013',...
'2014-04-10-2/data000',...
'2012-09-06-0/streamed/data001/data001',...
'2013-08-19-4/streamed/data001/data001',...
'2013-10-10-1/streamed/data000/data000'
};

cone_folders={...
    'streamed_data007-localmax-t9',...
    'data002-bayes-msf_5.00-BW-2-4',...
    'streamed_data001-localmax-t16',...
    'streamed_data013-bayes-msf_15.00-BW-2-4',...
    'streamed_data000-bayes-msf_10.00-all_BW-2-8',...
    'streamed_data001-localmax-t13',...
    'streamed_data001-bayes-msf_20.00-BW-3-4',...
    'streamed_data000-bayes-msf_10.00-BW-1-8'
    };

single_cone_eccentricity=[15,11,8.1,7.9,1,9,15.1,8.2];

for i_dataset=1:length(single_cone_datasets)
    
    disp(single_cone_datasets{i_dataset})
    
    % Find cone data
    datarun=load_data(single_cone_datasets{i_dataset});
    datarun=load_neurons(datarun);
    datarun=load_params(datarun);
    datarun=load_cones_njb(datarun, cone_folders{i_dataset});
    
    % find on and off parasols
    try
        [indices{1},~,~ ]=get_cell_indices(datarun, 'ON parasol');
    catch 
        warning('No ON parasols')
    end
    try
        [indices{2},~,~ ]=get_cell_indices(datarun, 'OFF parasol');
    catch
        warning('No OFF parasols')
    end
    
    ON_cones{i_dataset}=zeros(length(indices{1}),1);
    OFF_cones{i_dataset}=zeros(length(indices{2}),1);
    
    % find # of cones for each
    for i_cell=1:length(indices{1})
        cone_weights=datarun.cones.weights(:,indices{1}(i_cell));
        ON_cones{i_dataset}(i_cell)=sum(cone_weights>max(cone_weights)*thresh);
        clear cone_weights
    end
    for i_cell=1:length(indices{2})
        cone_weights=datarun.cones.weights(:,indices{2}(i_cell));
        OFF_cones{i_dataset}(i_cell)=sum(cone_weights>max(cone_weights)*thresh);
        clear cone_weights
    end
    
    clear indices
        
end

if 0
% save in some data structure
cone.ON.mean(i_dataset)=mean(ON_cones);
cone.ON.stderror(i_dataset)=std(ON_cones)/sqrt(length(ON_cones));
cone.OFF.mean(i_dataset)=mean(OFF_cones);
cone.OFF.stderror(i_dataset)=std(OFF_cones)/sqrt(length(OFF_cones));

% errorbar(single_cone_eccentricity, cone.ON.mean, cone.ON.stderror,'b')
% hold on
plot(single_cone_eccentricity, cone.OFF.mean,'ok')
ylabel('Average number of cones')
xlabel('Eccentricity (mm)')
% legend('On','Off')
end

%%
hold on
Colorset=varycolor(10);
for i_dataset=1:length(single_cone_datasets)
    % plot(single_cone_eccentricity(i_dataset)*ones(length(ON_cones{i_dataset}),1),ON_cones{i_dataset}, '.', 'Color', Colorset(i_dataset,:));
    plot(single_cone_eccentricity(i_dataset)*ones(length(OFF_cones{i_dataset}),1),OFF_cones{i_dataset}, '.', 'Color', Colorset(i_dataset,:));
end

