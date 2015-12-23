%plots pdf overview with rf-summary-classification
%greschner

file_path='/marte/snle/lab/Experiments/Array/Analysis/dataset_overview/';
if 1
    file_name='dataset_overview.txt';
    ext_name='';
else
    file_name='dataset_overview_preliminary.txt';
    ext_name='_preliminary';    
end

% open and load file
fid=fopen([file_path file_name]);
file=textscan(fid, '%s');
fclose(fid);
pathfile=file{1};

for i=22:length(pathfile)
    display(sprintf('%d - %d',i,length(pathfile)))
    
    clear datarun
    t=dir(['/' pathfile{i} '/d*.params']);
    if isempty(t)
        warning(['/' pathfile{i} '/*.params not found'])
    else
        datarun.names.rrs_params_path=['/' pathfile{i} '/' t.name];
        t=dir(['/' pathfile{i} '/rf-summary-classification.txt']);
        if ~isempty(t)
            datarun.names.rrs_classification_path=['/' pathfile{i} '/' t.name];
        end

        opt=struct('verbose',1,'load_params',1,'load_neurons',1);
        datarun=load_data(datarun,opt);     
        datarun=load_vision_cell_types(datarun);
        info (datarun)

        paramsFile = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun.names.rrs_params_path);

        fig=figure('Visible','off','PaperPosition',[0 0 12 8],'PaperSize',[12 8])
        set(gcf,'color','white');

        class_list=[1:4];
        for ii=5:length(datarun.cell_types)
            if ~isempty(datarun.cell_types{ii}.name)
                if isempty(strfind(datarun.cell_types{ii}.name,'duplicate'))&isempty(strfind(datarun.cell_types{ii}.name,'crap'))&isempty(strfind(datarun.cell_types{ii}.name,'weak'))&isempty(strfind(datarun.cell_types{ii}.name,'contaminated'))&isempty(strfind(datarun.cell_types{ii}.name,'edge'))    
                     class_list=[class_list ii];
                end
            end
        end
        
        ind=get_cell_indices(datarun,num2cell(class_list));
        t=zeros(length(ind),2);
        for ii=1:length(ind)
            t(ii,:)=datarun.vision.sta_fits{ind(ii)}.mean;
        end
        tt=sort(t(:,1));
        x_range=[tt(2)-(tt(end-1)-tt(2))*.1 tt(end-1)+(tt(end-1)+tt(2))*.1]; 
        %x_range=[min(t(:,1)) max(t(:,1))];
        tt=sort(t(:,2));
        y_range=[tt(2)-(tt(end-1)-tt(2))*.1 tt(end-1)+(tt(end-1)+tt(2))*.1];
        %y_range=[min(t(:,2)) max(t(:,2))];

        for ii=1:length(class_list)

            index=get_cell_indices(datarun,{class_list(ii)});
            cell_id=datarun.cell_ids(index);

            if ~isempty(cell_id)
                ttime=zeros(length(cell_id),3,length(paramsFile.getArrayCell(cell_id(1),'GreenTimeCourse')));
                tauto=zeros(length(cell_id),length(paramsFile.getArrayCell(cell_id(1),'Auto')));
                for iii=1:length(cell_id)
                    t=paramsFile.getArrayCell(cell_id(iii), 'RedTimeCourse');
                    if ~isempty(t)
                        ttime(iii,1,:)=t;
                        ttime(iii,2,:)=paramsFile.getArrayCell(cell_id(iii), 'GreenTimeCourse');
                        ttime(iii,3,:)=paramsFile.getArrayCell(cell_id(iii), 'BlueTimeCourse');
                        tauto(iii,:)=paramsFile.getArrayCell(cell_id(iii), 'Auto');
                    end
                end
                if length(cell_id)>1
                    for iiii=1:3
                        t=sum(abs(mean(ttime(:,iiii,:))));
                        for iii=1:length(cell_id)
                            ttime(iii,iiii,:)=ttime(iii,iiii,:)/sum(abs(ttime(iii,iiii,:)))*t;
                        end
                    end
                    t=sum(abs(mean(tauto)));
                    for iii=1:length(cell_id)
                        tauto(iii,:)=tauto(iii,:)/sum(abs(tauto(iii,:)))*t;
                    end
                end
               

                subplot(ceil(length(class_list)/2),6,(ii-1)*3+1)
                    plot_rf_fit(datarun, {class_list(ii)},'fill',false,'edge_color',[0 0 0]);
                    set(gca,'DataAspectRatio',[1 1 1]);
                    %set(gca,'PlotBoxAspectRatio',[1 1 1]);    
                    %set(gca,'XTick',[],'YTick',[],'Visible','off');
                    set(gca,'XTick',[],'YTick',[]);
                    set(gca,'xLim',x_range,'yLim',y_range);
                    title(datarun.cell_types{class_list(ii)}.name);
                    %ylabel(datarun.cell_types{class_list(ii)}.name);
                    box on
                subplot(ceil(length(class_list)/2),6,(ii-1)*3+2)
                    hold on
                    plot(squeeze(ttime(:,1,:))','r')
                    plot(squeeze(ttime(:,3,:))','b')
                    plot(squeeze(ttime(:,2,:))','g')
                    set(gca,'XTick',[],'YTick',[]);
                    set(gca,'xLim',[0  size(ttime,3)]);
                    box on;
                subplot(ceil(length(class_list)/2),6,(ii-1)*3+3)
                    plot(tauto','k')
                    set(gca,'XTick',[],'YTick',[]);
                    set(gca,'xLim',[0  size(tauto,2)]);
                    box on;
            end 
        end

        t=datarun.names.short_name;
        t=strrep(t,'_marte','');
        t=strrep(t,'_jacob','');
        t=strrep(t,'__braid_snle_analysis-archive_Experiments_Array_Analysis_','');
        if t(1)=='_', t=t(2:end); end
        print(fig,'-dpdf',sprintf('%s%s%s.pdf',file_path,t,ext_name));  
        close(fig)
    end
end

    
   
    
    
if 0
    fig = figure('Visible','off',...
      'PaperPosition',[0 0 6 4],...
      'PaperSize',[6 4])
  plot(rand(5))
  set(gca,'Position',[0 0 1 1])
  print(fig,'-dpdf',sprintf('%s%s.pdf',file_path,datarun.names.short_name))  
    
    
    % first, create a multipage ps
         fnam='foo.ps';
         delete(fnam);
    for i=1:10
         plot(rand(10,1));
         print('-dpsc2','-append',fnam);
    end
end