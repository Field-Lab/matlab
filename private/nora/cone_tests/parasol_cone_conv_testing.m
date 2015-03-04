thresh=0.3;
% 
% single_cone_datasets={...
% '2013-08-19-4'
% '2013-10-10-1',...
% '2014-04-10-2',...
% '2014-06-04-7'};

% '2012-09-13-2',... 5.5 mm

single_cone_datasets={...
'2011-05-11-6/data002'};

cone_folders={'data002-bayes-msf_5.00-BW-2-4'};

single_cone_eccentricity=[];

for i_dataset=1:1 %length(single_cone_datasets)
    
    datarun=load_data(single_cone_datasets{i_dataset});
    datarun=load_neurons(datarun);
    datarun=load_params(datarun);
    datarun=load_sta(datarun);
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

    ON_cones=zeros(length(indices{1}),1);
    OFF_cones=zeros(length(indices{2}),1);
    
    %% find # of cones for each
    for i_cell=1:length(indices{2})
        cone_weights=datarun.cones.weights(:,indices{2}(i_cell));
        sig_cones=cone_weights>max(cone_weights)*thresh;
        imagesc(datarun.stas.stas{indices{2}(i_cell)}(:,:,1,5)); colormap gray;
        hold on;
        plot(datarun.cones.centers(sig_cones,1),datarun.cones.centers(sig_cones,2),'+');
        hold off;
        pause()
        OFF_cones(i_cell)=sum(sig_cones);
        clear cone_weights
    end
    
    %%
    for i_cell=1:length(indices{2})
        cone_weights=datarun.cones.weights(:,indices{2}(i_cell));
        OFF_cones(i_cell)=sum(cone_weights>max(cone_weights)*thresh);
        clear cone_weights
    end
    
    % save in some data structure
    cone.ON.mean(i_dataset)=mean(ON_cones);
    cone.ON.stderror(i_dataset)=std(ON_cones)/sqrt(length(ON_cones));
    cone.OFF.mean(i_dataset)=mean(OFF_cones);
    cone.OFF.stderror(i_dataset)=std(OFF_cones)/sqrt(length(OFF_cones));

end