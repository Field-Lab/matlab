% Load Data
if 1
    % pick one of the datasets to load
    clear datarun
    
    if 0 %P23H data: age 525 days, no cells with RF
        datarun{1}.names.rrs_params_path='/Volumes/Brokedown/Analysis/Lauren/2009-10-23/data002/data002.params';
        datarun{1}.names.rrs_neurons_path='/Volumes/Brokedown/Analysis/Lauren/2009-10-23/data002/data002.neurons';
        datarun{1}.names.rrs_ei_path='/Volumes/Brokedown/Analysis/Lauren/2009-10-23/data002/data002.ei';
        datarun{1}.names.rrs_globals_path = '/Volumes/Brokedown/Analysis/Lauren/2009-10-23/data002/data002.globals';
        datarun{1}.out=[]; %specifies set of cells to exclude from figure
    end
    if 0 %P23H data: age 322 days, no ~40% of cells have RF
        datarun{1}.names.rrs_params_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data003/data003.params';
        datarun{1}.names.rrs_neurons_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data003/data003.neurons';
        datarun{1}.names.rrs_ei_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data003/data003.ei';
        datarun{1}.out=[]; %specifies set of cells to exclude from figure
    end
    if 0 %P23H data: age 322 days, no ~40% of cells have RF
        datarun{1}.names.rrs_params_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data004/data004.params';
        datarun{1}.names.rrs_neurons_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data004/data004.neurons';
        datarun{1}.names.rrs_ei_path='/Volumes/Brokedown/Analysis/Lauren/2007-11-26/data004/data004.ei';
        datarun{1}.out=[]; %specifies set of cells to exclude from figure
    end
    if 1 %P23H data: young animal used in first grant
        datarun{1}.names.rrs_params_path='/Volumes/Brokedown/Analysis/Lauren/2007-08-09/data004-cf-lh-complete/data004-cf-lh-complete.params';
        datarun{1}.names.rrs_neurons_path='/Volumes/Brokedown/Analysis/Lauren/2007-08-09/data004-cf-lh-complete/data004-cf-lh-complete.neurons';
        datarun{1}.names.rrs_ei_path='/Volumes/Brokedown/Analysis/Lauren/2007-08-09/data004-cf-lh-complete/data004-cf-lh-complete.ei';
        datarun{1}.names.rrs_globals_path = '/Volumes/Brokedown/Analysis/Lauren/2007-08-09/data004-cf-lh-complete/data004-cf-lh-complete.globals';
        datarun{1}.out=[]; %specifies set of cells to exclude from figure
        
        cellsInClasses{1} = [47 53 183 199 227 287 318 320 395 498 511 709 826 902]; %fast on transient
        cellsInClasses{2} = [62 167 182 256 316 406 436 496 541 616 646 707 721 766 796 830 932 946]; %fast off transient
        cellsInClasses{3} = [5 107 166 181 226 398 931]; %fast off sustained
        
    end
    
    
    if 0
        datarun{1}.names.rrs_params_path='/jacob/snle/lab/Experiments/Array/Analysis/2005-08-08-0/data010-4dot5/data010-4dot5.params';
        datarun{1}.names.rrs_neurons_path='/jacob/snle/lab/Experiments/Array/Analysis/2005-08-08-0/data010-4dot5/data010-4dot5.neurons';
        datarun{1}.names.rrs_ei_path='/jacob/snle/lab/Experiments/Array/Analysis/2005-08-08-0/data010-4dot5/data010-4dot5.ei';
        datarun{1}.out=[];
    end
    if 0
        datarun{1}.names.rrs_params_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-03-27-1/data008-gdf/data008.params';
        datarun{1}.names.rrs_neurons_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-03-27-1/data008-gdf/data008.neurons';
        datarun{1}.names.rrs_ei_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-03-27-1/data008-gdf/data008.ei';
        datarun{1}.out=[];
    end
    if 0
        datarun{1}.names.rrs_params_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-03-27-1/data008-gdf/data008.params';
        datarun{1}.names.rrs_neurons_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-03-27-1/data008-gdf/data008.neurons';
        datarun{1}.names.rrs_ei_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-03-27-1/data008-gdf/data008.ei';
        datarun{1}.out=[];
    end
    if 0
        datarun{1}.names.rrs_params_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-03-02-1/data016-gdf/data016/data016.params';
        datarun{1}.names.rrs_neurons_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-03-02-1/data016-gdf/data016/data016.neurons';
        datarun{1}.names.rrs_ei_path='/jacob/snle/lab/Experiments/Array/Analysis/2007-03-02-1/data016-gdf/data016/data016.ei';        
        datarun{1}.out=[];
    end

    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_sta',1,'load_ei',1);
    datarun=load_data(datarun,opt);
 
end 


%%

if 1 %load params    
    for j=1%:3
    
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
                datarun{j}.auto{i} = t/sqrt(sum(t.^2)); %normalizes acf by rms
                %datarun{j}.auto{i} = t;

                if 1 %use waveform from electrode with max ei amplitude
                    [junk t]=max(max(abs(datarun{j}.ei.eis{i}),[],2));
                    tt=datarun{j}.ei.eis{i}(t,:);
                else
                    [junk t]=sort(max(abs(datarun{j}.ei.eis{i}),[],2));
                    tt=mean(datarun{j}.ei.eis{i}(t(end-2:end),:));
                end
                %tt=tt/sum(abs(tt));
                tt=tt/sqrt(sum(tt.^2)); %normalizes ei waveform by rms
                %tt=tt/sum(abs(tt(find(tt<0))));
                [junk t]=min(tt);
                datarun{j}.ei_waveform{i}=tt(t-8:t+37)';  %aligns ei waveform minima               
        end
    end
end

% if 1 %plot test
%     clf
%     hold on
%     box on
%     col=colormap(lines);
%     j=1
%     %for i=[4 3 1 2]
%         %index=get_cell_indices(datarun{j},{i});
%         index=get_cell_indices(datarun{j},'all');
%         for ii=1:length(index)
%             datarun{j}.cell_ids(index(ii));
%             %plot(datarun{j}.ei_waveform{index(ii)},'color',col(i,:)) 
%             plot(datarun{j}.auto{index(ii)},'color',col(i,:)) 
%             %plot(datarun{j}.green_time_course{index(ii)}) 
%         end
%     %end
% end

%% make plot based on non-visual paramaters
if 1 %plot PCA   FIGURE GRANT
    datarun{1}.out=[];
    datarun{2}.out=[];
    datarun{3}.out=[];
    
    
    keyboard
    
    clf
    col=colormap(lines);
    for j=1
        r=[];rr=[];rrr=[];c=[];id=[];
        %for i=[1 2 3 4]
        for i = [1 2 3];
            %index=get_cell_indices(datarun{j},{i});
            index = get_cell_indices(datarun{j}, cellsInClasses{i});
            
            for ii=1:length(index)
                if ~ismember(datarun{j}.cell_ids(index(ii)), datarun{j}.out)
                    id=[id datarun{j}.cell_ids(index(ii))];
                    c=[c; col(i,:)];

                    %r=[r datarun{j}.green_time_course{index(ii)}];
                    
                    r=[r datarun{j}.auto{index(ii)}(4:8)]; %one window of acf
                    rr=[rr datarun{j}.auto{index(ii)}(8:end)]; %another window of acf
                    %rr=[rr datarun{j}.auto{index(ii)}];

                    rrr=[rrr datarun{j}.ei_waveform{index(ii)}(10:end)]; %ei waveform
                  
                end
            end
        end
        [coeff,scores1] = princomp(r');
        [coeff,scores2] = princomp(rr');
        [coeff,scores3] = princomp(rrr');
        
        %what do these do?
        %scores4=mass(rr(1:end,:))';
        scores5=sum(rr(1:10,:))'./sum(rr(20:30,:))';
   
        %subplot(3,2,(j-1)*2+2)
            
            for ii=1:length(scores2(:,1))
                if ~ismember(id(ii),[6831 1817])%6831 
                    plot(scores1(ii,1),scores2(ii,1),'.','color',c(ii,:),'MarkerSize',12)%12
                    %plot(scores5(ii,1),scores4(ii,1),'.','color',c(ii,:),'MarkerSize',5)%12
                    %plot(scores2(ii,1),scores3(ii,1),'.','color',c(ii,:),'MarkerSize',5)%12
                    %plot3(scores2(ii,1),scores2(ii,2),scores3(ii,1),'.','color',c(ii,:),'MarkerSize',5)%12
                    hold on
                end
            end
            set(gca,'PlotBoxAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);
            %set(gca,'PlotBoxAspectRatio',[10 16 1],'Xtick',[],'Ytick',[]);
            %[xrange yrange]=autoscale(scores3(:,1),scores2(:,1),'border', .15);
            %set(gca,'XLim',xrange,'YLim',yrange);
            box on

    end
end

%% make plot based on visual parameters

if 0 %plot PCA   FIGURE GRANT
    datarun{1}.out=[];
    datarun{2}.out=[];
    datarun{3}.out=[];
    
    
    clf
    col=colormap(lines);
    for j=1
        r=[];rr=[];rrr=[];c=[];id=[];
        for i=[1 2 3 4]
            index=get_cell_indices(datarun{j},{i});
            
            for ii=1:length(index)
                if ~ismember(datarun{j}.cell_ids(index(ii)), datarun{j}.out)
                    id=[id datarun{j}.cell_ids(index(ii))];
                    c=[c; col(i,:)];

                    r=[r ext(datarun{j}.green_time_course{index(ii)})];
                    
                    rr=[rr pi*prod(datarun{j}.vision.sta_fits{index(ii)}.sd)*8*5.75];
                   
                end
            end
        end
  
        %subplot(3,2,(j-1)*2+2)
            
            for ii=1:length(r)
                if ~ismember(id(ii),[6831])%6831 1817
                    plot(r(ii),rr(ii),'.','color',c(ii,:),'MarkerSize',12)%12
                    
                    hold on
                end
            end
            set(gca,'PlotBoxAspectRatio',[1 1 1],'Xtick',[],'Ytick',[]);
            %[xrange yrange]=autoscale(scores3(:,1),scores2(:,1),'border', .15);
            %set(gca,'XLim',xrange,'YLim',yrange);
            box on

    end
end
%[t tt]=min(abs(scores2(:,1)-(-.2368)));id(tt)




%LOAD

% save('/snle/home/greschner/Desktop/AC/ei_2007-03-02-1data016-gdf.mat', 'datarun')

if 0 %select electrode
datarun{1}=load_ei(datarun{1}, 'all','array_type',512)       

    params.cutoff=0.03; 
    params.field_name='marks'; 
    params.overwrite_name=0;

    nr=1;
    cell_specification={3};

    %get marks  
        if ~params.overwrite_name & isfield(datarun{nr}.ei, params.field_name);
            ei_marks=datarun{nr}.ei.(params.field_name);
        else
            ei_marks{length(datarun{nr}.cell_ids)}=[];
        end

    % get cell numbers
        index=get_cell_indices(datarun{nr},cell_specification);


        disp('When finished positioning polygon double-clicking or right-clicking inside the polygon')

        clf
        r=[min(datarun{nr}.ei.position); max(datarun{nr}.ei.position)]*1.1;
        set(gca,'XLim',r(:,1),'YLim',r(:,2),'YTick',[],'XTick',[],'PlotBoxAspectRatio',[2 1 1]);
        box on;

        k=1; 
        kmin=1; 
        kmax=length(index);
        if kmax>1
            handel1= uicontrol(gcf,...
            'Style','slider',...
            'Min' ,kmin,'Max',kmax, ...
            'Units','normalized', ...
            'Position',[0,0,.6,.04], ...
            'Value', k,...
            'SliderStep',[1/(kmax-kmin) 1/(kmax-kmin)],...
            'CallBack', 'uiresume;');
        end
        handel2= uicontrol(gcf,...
            'Style','pushbutton', ...
            'Units','normalized', ...
            'Position',[0.61 0 .06 .04], ...
            'String','del',...
            'Callback','marks=[]; hold off; image(0); hold on; plot_ei(datarun{nr},datarun{nr}.cell_ids(index(k))); uiresume;');


        while k

            k=round(get(handel1,'Value'));

            image(0);
            hold on;


            %plot_ei(datarun{nr},datarun{nr}.cell_ids(index(k)),0,'cutoff',params.cutoff,'neg_color',[0 0 0],'pos_color',[0 0 0],'fliplr',1,'flipud',1,'highlight',datarun{nr}.channels(index(k)));
            plot_ei(datarun{nr},datarun{nr}.cell_ids(index(k)),'cutoff',params.cutoff,'neg_color',[0 0 0],'pos_color',[0 0 0],'fliplr',1,'flipud',1,'highlight',datarun{nr}.channels(index(k)));
 
            title(sprintf('cell-id=%d %d/%d',datarun{nr}.cell_ids(index(k)),k,length(index))); 

            marks=ei_marks{index(k)};
            if ~isempty(marks)
                plot(datarun{nr}.ei.position(marks,1), datarun{nr}.ei.position(marks,2),'r.');
            end    

            [junk,xi,yi]=roipoly;
            t1=inpolygon(datarun{nr}.ei.position(:,1),datarun{nr}.ei.position(:,2),xi,yi);

            ei=datarun{nr}.ei.eis{index(k)};
            ei=ei/max(abs(ei(:)));
            t2=max(abs(ei),[],2);

            marks=[marks; find(t1 & t2>=params.cutoff)];
            plot(datarun{nr}.ei.position(marks,1), datarun{nr}.ei.position(marks,2),'r.');

            ei_marks{index(k)}=marks;
            hold off

            %assignin('base','ei_marks',ei_marks);
            datarun{nr}.ei.(params.field_name)=ei_marks;

            uiwait; 
        end

  ei_marks=datarun{1}.ei.marks; save('/snle/home/greschner/Desktop/AC/ei_2009-12-03-2-003-nwpca.mat', 'ei_marks');

end    
    
if 0
    p=datarun{1}.ei.position;

    x=min(p(:,1)):2:max(p(:,1));
    y=min(p(:,2)):2:max(p(:,2));

    position_fine=zeros(length(x),length(y));
    for i=1:length(x)
        for ii=1:length(y)
            t=(p(:,1)-x(i)).^2+(p(:,2)-y(ii)).^2;
            [junk,position_fine(i,ii)]=min(t);
        end
    end

    plot_mat(position_fine(1:200,1:200))
end



if 0 %get ei mosaic cont
    f=fspecial('gaussian',100,15);
    thr1=.35;
    
    index=get_cell_indices(datarun{1},{3});   
    for i=1:length(index)

        [junk t]=max(max(abs(datarun{1}.ei.eis{i}),[],2));
        tt=datarun{j}.ei.eis{i}(t,:);
        tnorm=sqrt(sum(tt.^2));
        
        m=datarun{1}.ei.marks{index(i)};
        if ~isempty(m)
            p=zeros(size(position_fine));
            for ii=1:length(m)
                p(find(position_fine==m(ii)))=max(abs(datarun{1}.ei.eis{index(i)}(m(ii),:)/tnorm));
            end

            p=p';
            p=imfilter(p,f,'replicate');

            t=thr1*max(p(:));
            h1=contourc(p,[t t]);        
            %get longest conture
            ii=1;
            n=1;
            k=ii;
            while ii<length(h1)
                if h1(2,ii)>n
                    n=h1(2,ii);
                    k=ii+1;
                end
                ii=ii+h1(2,ii)+1;
            end

            t=[h1(:,k:k+h1(2,k-1)-1),h1(:,k)];
            datarun{1}.ei.cont{index(i)}=t;

            if 0
                clf
                plot_mat(p)
                hold on
                plot(t(1,:),t(2,:),'k');
            end
        else
            missed=i
        end
    end
    
end
if 0 %plot ei mosaic
    clf  
    hold on

    %index=get_cell_indices(datarun{1},{1});   
    for i=1:length(index)
        
        c=datarun{1}.ei.cont{index(i)};  
        if ~isempty(c)
            plot(datarun{1}.ei.cont{index(i)}(1,:),datarun{1}.ei.cont{index(i)}(2,:),'k');
        end
   
    end
    set(gca,'YTick',[],'XTick',[],'YLim',[1,size(position_fine,2)],'XLim',[1,size(position_fine,1)],'PlotBoxAspectRatio',[2 1 1]);
    box on
end
if 0
    clf
    plot_rf_fit(datarun{1}, {4},'sd_radius',1.1);
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gca,'PlotBoxAspectRatio',[1 1 1]); 
    set(gca,'XTick',[],'YTick',[]);
    set(gca,'XLim',[0 80],'YLim',[0 40]);
    box on;
end















