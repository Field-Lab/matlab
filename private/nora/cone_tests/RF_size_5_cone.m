
single_cone_datasets={...
'2012-09-21-2/data007',...
'2011-05-11-6/data002',...
'2012-09-24-5/data001',...
'2011-06-24-6/data013',...
'2012-09-06-0/streamed/data001/data001',...
'2013-10-10-1/streamed/data000/data000'
};

cone_folders={...
    'streamed_data007-localmax-t9',...
    'data002-bayes-msf_5.00-BW-2-4',...
    'streamed_data001-localmax-t16',...
    'streamed_data013-bayes-msf_15.00-BW-2-4',...
    'streamed_data001-localmax-t13',...
    'streamed_data000-bayes-msf_10.00-BW-1-8'
    };

single_cone_eccentricity=[15,11,8.1,7.9,9,8.2];


for i_dataset=1:length(single_cone_datasets)
    
    count=1;
    
    datarun=load_data(single_cone_datasets{i_dataset});
    % datarun=load_neurons(datarun);
    datarun=load_params(datarun);
    % datarun=load_sta(datarun);
    datarun=load_cones_njb(datarun, cone_folders{i_dataset});
    n_cones=length(datarun.cones.types);
    
    plot(datarun.cones.centers(:,1), datarun.cones.centers(:,2), '.k')
    [x,y]=ginput();
    
    for i_click=1:length(x)
        mean_cones=datarun.cones.centers-repmat([x(i_click), y(i_click)],n_cones,1);
        norm_mean_cones=diag(mean_cones*mean_cones');
        temp=sort(norm_mean_cones);
        RF_radius{i_dataset}(count)=sqrt(temp(5))*5;
        ecc{i_dataset}(count)=single_cone_eccentricity(i_dataset);
        count=count+1;
    end
    
    clear x y
  
end

%%
figure;
hold on
Colorset=varycolor(10);
for i_dataset=1:length(single_cone_datasets)
    plot(single_cone_eccentricity(i_dataset)*ones(length(RF_radius{i_dataset}),1),RF_radius{i_dataset}, '.', 'Color', Colorset(i_dataset,:));
end