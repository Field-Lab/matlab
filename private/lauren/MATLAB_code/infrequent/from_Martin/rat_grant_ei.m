% Load Data
if 1
    clear datarun
    
    if 1
        datarun{1}.names.rrs_params_path='/snle/lab/Experiments/Array/Analysis/2009-09-03-1/data002/data002.params';
        datarun{1}.names.rrs_neurons_path='/snle/lab/Experiments/Array/Analysis/2009-09-03-1/data002/data002.neurons';
        datarun{1}.names.rrs_ei_path='/snle/lab/Experiments/Array/Analysis/2009-09-03-1/data002/data002.ei';
        datarun{1}.names.rrs_globals_path='/snle/lab/Experiments/Array/Analysis/2009-09-03-1/data002/data002.globals';
        %datarun{1}.celltypes=[13 21 11];
        datarun{1}.celltypes = [14 21 12];
        datarun{1}.out=[829  663 766 931  586 706 811];
    end
    if 1
        datarun{2}.names.rrs_params_path='/snle/lab/Experiments/Array/Analysis/2009-09-06-0/data001/data001.params';
        datarun{2}.names.rrs_neurons_path='/snle/lab/Experiments/Array/Analysis/2009-09-06-0/data001/data001.neurons';
        datarun{2}.names.rrs_ei_path='/snle/lab/Experiments/Array/Analysis/2009-09-06-0/data001/data001.ei';
        datarun{2}.names.rrs_globals_path='/snle/lab/Experiments/Array/Analysis/2009-09-06-0/data001/data001.globals';
        %datarun{2}.celltypes=[8 19 14];
        datarun{2}.celltypes = [9 20 15];
        datarun{2}.out=[436 586 811  136 601 706  64 901];
    end
    if 1
        datarun{3}.names.rrs_params_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-08-09/data004-cf-lh-complete/data004-cf-lh-complete.params';
        datarun{3}.names.rrs_neurons_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-08-09/data004-cf-lh-complete/data004-cf-lh-complete.neurons';
        datarun{3}.names.rrs_ei_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-08-09/data004-cf-lh-complete/data004-cf-lh-complete.ei';
        datarun{3}.names.rrs_globals_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-08-09/data004-cf-lh-complete/data004-cf-lh-complete.globals';
        %datarun{3}.celltypes=[13 12 14];
        datarun{3}.celltypes=[14 13 15];
        datarun{3}.out=[5  167 182 616 646 707 766  47 287];
    end
    
    if 0 %P23H data: age 525 days, no cells with RF
        datarun{4}.names.rrs_params_path='/Volumes/Brokedown/Analysis/Lauren/2009-10-23/data002/data002.params';
        datarun{4}.names.rrs_neurons_path='/Volumes/Brokedown/Analysis/Lauren/2009-10-23/data002/data002.neurons';
        datarun{4}.names.rrs_ei_path='/Volumes/Brokedown/Analysis/Lauren/2009-10-23/data002/data002.ei';
        datarun{4}.names.rrs_globals_path = '/Volumes/Brokedown/Analysis/Lauren/2009-10-23/data002/data002.globals';
        datarun{4}.celltypes=[13];
        datarun{4}.out=[]; %specifies set of cells to exclude from figure
    end
    if 0 %P23H data: age 322 days, ~40% of cells have RF
        datarun{4}.names.rrs_params_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data003/data003.params';
        datarun{4}.names.rrs_neurons_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data003/data003.neurons';
        datarun{4}.names.rrs_ei_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data003/data003.ei';
        datarun{4}.names.rrs_globals_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data003/data003.globals';
        datarun{4}.celltypes=[8];
        datarun{4}.out=[]; %specifies set of cells to exclude from figure
    end
    if 0 %P23H data: age 322 days, ~40% of cells have RF
        datarun{4}.names.rrs_params_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data004/data004.params';
        datarun{4}.names.rrs_neurons_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data004/data004.neurons';
        datarun{4}.names.rrs_ei_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data004/data004.ei';
        datarun{4}.out=[]; %specifies set of cells to exclude from figure
    end
    
    
    

    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_sta',1,'load_ei',1);
    datarun=load_data(datarun,opt);

    %datarun{1}=load_ei(datarun{1}, 'all','array_type',61);  
    %datarun{2}=load_ei(datarun{2}, 'all','array_type',61);  
    %datarun{3}=load_ei(datarun{3}, 'all','array_type',61);  
end 

%datarun{j}.cell_types{datarun{j}.celltypes(i)}.name)


%%



if 0 %plot spike waveforms
    clf
    
    for j=1:3
        for i=1:3
            index=get_cell_indices(datarun{j},{datarun{j}.celltypes(i)});   

            r=zeros(length(index),size(datarun{j}.ei.eis{index(1)},2));
            for ii=1:length(index)
                [junk t]=max(max(abs(datarun{j}.ei.eis{index(ii)}),[],2));
                tt=datarun{j}.ei.eis{index(ii)}(t,:);
                r(ii,:)=tt/sum(abs(tt));
            end
            subplot(4,3,(i-1)*3+j)
                plot(r','b')
                hold on
                plot(mean(r),'r')
                set(gca,'XLim',[10 50]);
            subplot(4,3,3*3+j)
                plot(mean(r),'r')
                hold on
                set(gca,'XLim',[10 50]);
        end
    end
    
end

%%


if 0 %plot spike waveforms
    clf
    
    for j=1:3
        i=1;
            index=get_cell_indices(datarun{j},{datarun{j}.celltypes(i)});   

            r=zeros(length(index),size(datarun{j}.ei.eis{index(1)},2));
            for ii=1:length(index)
                [junk t]=max(max(abs(datarun{j}.ei.eis{index(ii)}),[],2));
                tt=datarun{j}.ei.eis{index(ii)}(t,:);
                r(ii,:)=tt/sum(abs(tt));
            end
            subplot(3,1,j)
                plot(r','b')
                hold on
                plot(mean(r),'r.')
                set(gca,'XLim',[10 50]);
            
        i=2;
            index=get_cell_indices(datarun{j},{datarun{j}.celltypes(i)});   

            r=zeros(length(index),size(datarun{j}.ei.eis{index(1)},2));
            for ii=1:length(index)
                [junk t]=max(max(abs(datarun{j}.ei.eis{index(ii)}),[],2));
                tt=datarun{j}.ei.eis{index(ii)}(t,:);
                r(ii,:)=tt/sum(abs(tt));
            end
            subplot(3,1,j)
                plot(r','g')
                hold on
                plot(mean(r),'r*')
                set(gca,'XLim',[10 50]);
    end
    
end


%%


if 0 %load ei
    clear res
    %for j=1:3
    for j=1:4
        %for i=1:3
        for i = 1:length(datarun{j}.celltypes)
            index=get_cell_indices(datarun{j},{datarun{j}.celltypes(i)});   
            

            r=zeros(length(index),50);
            for ii=1:length(index)
                [junk t]=max(max(abs(datarun{j}.ei.eis{index(ii)}),[],2));
                tt=datarun{j}.ei.eis{index(ii)}(t,:);
                tt=tt/sum(abs(tt));
                %tt=tt/sum(abs(tt(find(tt<0))));
                
                [junk t]=min(tt);
                r(ii,:)=tt(t-10:t+39);
            end
            
            res{j,i}.r=r;
        end
    end   
end

if 0 %plot ei
    clf
    col=colormap(lines);
    
    for j=1:3
        for i=1:3
            subplot(3,4,(j-1)*4+i)
                plot(res{j,i}.r','color',col(i,:))
                hold on
               % plot(mean(res{j,i}.r),'k.-')
                set(gca,'XLim',[1 50]);
                title(sprintf('%d - %s',j,datarun{j}.cell_types{datarun{j}.celltypes(i)}.name))
            subplot(3,4,(j-1)*4+4)
                hold on
                plot(mean(res{j,i}.r),'color',col(i,:))
                set(gca,'XLim',[1 50]);
        end
    end
    
end

if 0 %plot pca from ei scores
    j=3;
    t=[res{j,1}.r; res{j,2}.r; res{j,3}.r];
    [coeff,scores] = princomp(t);
    clf
    hold on
    plot( scores(1:size(res{j,1}.r,1),1) , scores(1:size(res{j,1}.r,1),2) ,'.')
    plot( scores(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),1) , scores(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),2),'g.')
    plot( scores(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:end,1) , scores(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:end,2),'r.')
end

%%


if 0 %load params
    for j=1:3
    
        paramsFile = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun{j}.names.rrs_params_path);

        for i=1:length(datarun{j}.cell_ids)

            t = paramsFile.getArrayCell(datarun{j}.cell_ids(i), 'GreenTimeCourse');
                datarun{j}.green_time_course{i} = t/sum(abs(t));
            t = paramsFile.getArrayCell(datarun{j}.cell_ids(i), 'Auto');
                datarun{j}.auto{i} = t/sum(abs(t));

        end
    end
end




if 1 %load params -- this version used in grant submission
    for j=1:3
    %for j = 1:4
    
        paramsFile = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun{j}.names.rrs_params_path);

        for i=1:length(datarun{j}.cell_ids)

            t = paramsFile.getArrayCell(datarun{j}.cell_ids(i), 'GreenTimeCourse');
            t1 = paramsFile.getArrayCell(datarun{j}.cell_ids(i), 'RedTimeCourse');
            t2 = paramsFile.getArrayCell(datarun{j}.cell_ids(i), 'BlueTimeCourse');
                %datarun{j}.green_time_course{i} = t/sum(abs(t));
                %datarun{j}.green_time_course{i} = t/sum(abs(t+t1+t2));
                datarun{j}.green_time_course{i} = t/sqrt(sum(t.^2));
                %datarun{j}.green_time_course{i} = t/sqrt(sum((t+t1+t2).^2));
            t = paramsFile.getArrayCell(datarun{j}.cell_ids(i), 'Auto');
                %datarun{j}.auto{i} = t/sum(abs(t));
                datarun{j}.auto{i} = t/sqrt(sum(t.^2)); %used in grant
                %datarun{j}.auto{i} = t;

                [junk t]=max(max(abs(datarun{j}.ei.eis{i}),[],2));
                tt=datarun{j}.ei.eis{i}(t,:);
                tt=tt/sum(abs(tt));
                %tt=tt/sum(abs(tt(find(tt<0))));
                [junk t]=min(tt);
                datarun{j}.ei_waveform{i} =tt(t-10:t+39);                 
        end
    end
end
if 0 %plot test
    clf
    hold on
    j=1
    for i=[3]
        index=get_cell_indices(datarun{j},{datarun{j}.celltypes(i)});
        for ii=1:length(index)
            datarun{j}.cell_ids(index(ii))
            %plot(datarun{j}.auto{index(ii)}) 
            plot(datarun{j}.green_time_course{index(ii)}) 
        end
    end
end


if 0 %
    clear res
    for j=1:3
        for i=1:3
            index=get_cell_indices(datarun{j},{datarun{j}.celltypes(i)});   

            r=zeros(length(index),50);
            for ii=1:length(index)
                [junk t]=max(max(abs(datarun{j}.ei.eis{index(ii)}),[],2));
                tt=datarun{j}.ei.eis{index(ii)}(t,:);
                tt=tt/sum(abs(tt));
                %tt=tt/sum(abs(tt(find(tt<0))));
                
                [junk t]=min(tt);
                r(ii,:)=tt(t-10:t+39);
            end
            
            res{j,i}.r=r;
        end
    end   
end


if 0 %plot PCA   
    clf
    for j=1:3;
        r=[];rr=[];c=[];
        for i=1:3
            index=get_cell_indices(datarun{j},{datarun{j}.celltypes(i)});   
            for ii=1:length(index)
                c=[c datarun{j}.cell_ids(index(ii))];
                r=[r datarun{j}.auto{index(ii)}(70:end)];
                rr=[rr datarun{j}.green_time_course{index(ii)}];
            end
        end
        [coeff,scores1] = princomp(rr');
        [coeff,scores2] = princomp(r');

        %hack
        if j>1
            scores1=-scores1;
        end
        
        subplot(3,2,(j-1)*2+1)
        hold on
        col=colormap(lines);
        plot( scores1(1:size(res{j,1}.r,1),1) , scores2(1:size(res{j,1}.r,1),1) ,'.','color',col(1,:),'MarkerSize',12)
        plot( scores1(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),1) , scores2(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),1),'.','color',col(2,:),'MarkerSize',12)
        plot( scores1(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:end,1) , scores2(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:end,1),'.','color',col(3,:),'MarkerSize',12)
        set(gca,'PlotBoxAspectRatio',[10 16 1],'Xtick',[],'Ytick',[]);
        box on
    end
    
    for j=1:3;
        t=[res{j,1}.r; res{j,2}.r; res{j,3}.r];
        [coeff,scores1] = princomp(t);

        subplot(3,2,(j-1)*2+2)
        hold on
        col=colormap(lines);
        plot( scores1(1:size(res{j,1}.r,1),1) , scores2(1:size(res{j,1}.r,1),1) ,'.','color',col(1,:),'MarkerSize',12)
        plot( scores1(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),1) , scores2(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),1),'.','color',col(2,:),'MarkerSize',12)
        plot( scores1(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:end,1) , scores2(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:end,1),'.','color',col(3,:),'MarkerSize',12)
        set(gca,'PlotBoxAspectRatio',[10 16 1],'Xtick',[],'Ytick',[]);
        box on
    end
end

%%

if 1 %plot PCA   FIGURE GRANT (used in grant submission)
    datarun{1}.out=[];
    datarun{2}.out=[];
    datarun{3}.out=[395 53]; %53 has contaminated by an axon, 395 on edge of array
    datarun{4}.out = [];
    
    
    clf
    col=colormap(lines);
    for j=1:3;
    %for j = 1:4 %loops through datasets
        r=[];rr=[];rrr=[];c=[];id=[];
        %for i=1:3 %loops through cell types
        for i = 1:length(datarun{j}.celltypes)
            index=get_cell_indices(datarun{j},{datarun{j}.celltypes(i)});   
            for ii=1:length(index)
                if ~ismember(datarun{j}.cell_ids(index(ii)), datarun{j}.out)
                    id=[id datarun{j}.cell_ids(index(ii))];
                    c=[c; col(i,:)];

                    r=[r datarun{j}.green_time_course{index(ii)}];
                    %rr=[rr datarun{j}.auto{index(ii)}(70:end)];
                    rr=[rr datarun{j}.auto{index(ii)}(20:end)];
                    %rr=[rr datarun{j}.auto{index(ii)}(5:end)]; %works well for late-stage P23H
                    %rr=[rr datarun{j}.auto{index(ii)}(12:end)]; %good compromise for all pieces
                    

                    [junk t]=max(max(abs(datarun{j}.ei.eis{index(ii)}),[],2));
                    tt=datarun{j}.ei.eis{index(ii)}(t,:);
                    tt=tt/sum(abs(tt));
                    %tt=tt/sqrt(sum(tt.^2));
                    [junk t]=min(tt);
                    rrr=[rrr tt(t-10:t+39)];
                end
            end
        end
        [coeff,scores1] = princomp(r');
        [coeff,scores2] = princomp(rr');
        [coeff,scores3] = princomp(rrr');
        
        scores2save{j} = scores2;

        %hack
        if j>1
            scores1=-scores1;
        end
        
        %subplot(3,2,(j-1)*2+1)
        subplot(4,2,(j-1)*2+1)
            hold on
            for ii=1:length(scores1(:,1))
                plot(scores1(ii,1),scores2(ii,1),'.','color',c(ii,:),'MarkerSize',12)
            end
            set(gca,'PlotBoxAspectRatio',[10 16 1],'Xtick',[],'Ytick',[]);
            %set(gca,'PlotBoxAspectRatio',[10 10 1]);
            [xrange yrange]=autoscale(scores1(:,1),scores2(:,1),'border', .15);
            set(gca,'XLim',xrange,'YLim',yrange);
            box on

    end
    
    % out
    %datarun{1}.out=[829  663 766 931  586 706 811];
    %datarun{2}.out=[436 586 811  136 601 706  64 901];
    %datarun{3}.out=[5  167 182 616 646 707 766  47 287];
    
    %datarun{1}.out=[];
    %datarun{2}.out=[];
    %datarun{3}.out=[395];
    

    for j=1:3;
    %for j = 1:4
        r=[];rr=[];rrr=[];c=[];id=[];
        %for i=1:3
        for i = 1:length(datarun{j}.celltypes)
            index=get_cell_indices(datarun{j},{datarun{j}.celltypes(i)});
            
            for ii=1:length(index)
                if ~ismember(datarun{j}.cell_ids(index(ii)), datarun{j}.out)
                    id=[id datarun{j}.cell_ids(index(ii))];
                    c=[c; col(i,:)];

                    r=[r datarun{j}.green_time_course{index(ii)}];
                    %rr=[rr datarun{j}.auto{index(ii)}(70:end)];
                    %rr=[rr datarun{j}.auto{index(ii)}];
                    %rr=[rr datarun{j}.auto{index(ii)}(20:end)];

                    [junk t]=max(max(abs(datarun{j}.ei.eis{index(ii)}),[],2));
                    tt=datarun{j}.ei.eis{index(ii)}(t,:);
                    %tt=tt/sum(abs(tt));
                    tt=tt/sqrt(sum(tt.^2));
                    [junk t]=min(tt);
                    rrr=[rrr tt(t-10+0:t+39)'];
                end
            end
        end
        %[coeff,scores2] = princomp(rr'); %acf
        scores2 = scores2save{j};
        [coeff,scores3] = princomp(rrr'); %spike waveform
   
        %subplot(3,2,(j-1)*2+2)
        subplot(4,2,(j-1)*2+2)
            hold on
            for ii=1:length(scores2(:,1))
                %plot(scores3(ii,1),scores2(ii,1),'.','color',c(ii,:),'MarkerSize',12)
                plot(scores2(ii,2),scores2(ii,1),'.','color',c(ii,:),'MarkerSize',12)
            end
            set(gca,'PlotBoxAspectRatio',[10 16 1],'Xtick',[],'Ytick',[]);
            %set(gca,'PlotBoxAspectRatio',[10 10 1]);
            %[xrange yrange]=autoscale(scores3(:,1),scores2(:,1),'border', .15);
            [xrange yrange]=autoscale(scores2(:,2),scores2(:,1),'border', .15);
            set(gca,'XLim',xrange,'YLim',yrange);
            box on

    end
end


%%


if 0%plot3
    j=1;
    t=[res{j,1}.r; res{j,2}.r; res{j,3}.r];
    [coeff,scores1] = princomp(t);
    
    r=[];rr=[];
    for i=1:3
        index=get_cell_indices(datarun{j},{datarun{j}.celltypes(i)});   
        for ii=1:length(index)
            r=[r datarun{j}.auto{index(ii)}];
            rr=[rr datarun{j}.green_time_course{index(ii)}];
        end
    end
    [coeff,scores3] = princomp(rr');
    [coeff,scores2] = princomp(r');

    clf
    plot3( scores1(1:size(res{j,1}.r,1),1) , scores2(1:size(res{j,1}.r,1),1) , scores3(1:size(res{j,1}.r,1),1),'.')
    hold on
    plot3( scores1(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),1) , scores2(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),1), scores3(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),1),'g.')
    plot3( scores1(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:end,1) , scores2(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:end,1), scores3(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:end,1),'r.')

    %box on
end

if 0%pca agains all
    datarun{1}.celltypes_all=[13 24 11 8 9 12 14 15 16 17 19];
    datarun{2}.celltypes_all=[8 22 15 7 9 10 16 19 21 ];
    datarun{3}.celltypes_all=[13 12 14 8 10 17 18 19];
    
    j=1;
    t=[res{j,1}.r; res{j,2}.r; res{j,3}.r];
    
    
    r=[];rr=[];rrr=[];
    for i=1:length(datarun{j}.celltypes_all)
        index=get_cell_indices(datarun{j},{datarun{j}.celltypes_all(i)});   
        for ii=1:length(index)
            r=[r datarun{j}.auto{index(ii)}];
            rr=[rr datarun{j}.green_time_course{index(ii)}];
            rrr=[rrr datarun{j}.ei_waveform{index(ii)}'];
        end
    end
    [coeff,scores1] = princomp(rr');
    [coeff,scores2] = princomp(r');
    %[coeff,scores1] = princomp(rrr');

    clf
    hold on
    plot( scores1(1:size(res{j,1}.r,1),1) , scores2(1:size(res{j,1}.r,1),1) ,'.','MarkerSize',12)
    plot( scores1(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),1) , scores2(size(res{j,1}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1),1),'g.','MarkerSize',12)
    plot( scores1(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1)+size(res{j,3}.r,1),1) , scores2(size(res{j,1}.r,1)+size(res{j,2}.r,1)+1:size(res{j,1}.r,1)+size(res{j,2}.r,1)+size(res{j,3}.r,1),1),'r.','MarkerSize',12)
    plot( scores1(size(res{j,1}.r,1)+size(res{j,2}.r,1)+size(res{j,3}.r,1)+1:end,1) , scores2(size(res{j,1}.r,1)+size(res{j,2}.r,1)+size(res{j,3}.r,1)+1:end,1),'k.','MarkerSize',12)

    box on
end










