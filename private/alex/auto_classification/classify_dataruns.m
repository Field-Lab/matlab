datarun = load_data(fullfile(server_path(), '2015-01-29-0/data003/data003'));
datarun = load_params(datarun,'verbose',1);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_ei(datarun, 'all');
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);

if datarun.stimulus.independent=='t'
    rgb=1;
else
    rgb=0;
end


ncells=length(datarun.cell_ids);
t=[];pairs=[];cnt=1;p=[];cnt1=1;
for i=1:ncells-1
    for j=i+1:ncells
        t=corr2(datarun.ei.eis{i},datarun.ei.eis{j});
        p(cnt1)=t;cnt1=cnt1+1;
        if t>0.6
            pairs(cnt,:)=[i j];
            cnt=cnt+1;
        end
    end
end


dubl=datarun.cell_ids(pairs(:,1));
uncells=setdiff(datarun.cell_types{3}.cell_ids,dubl);

figure
plot_rf_summaries(datarun, uncells, 'clear', false,  'plot_fits', true, 'fit_color', 'r')



tc_length=length(datarun.vision.timecourses(1).r(16:end));
time_courses=zeros(tc_length,ncells);
isis=zeros(100,ncells);

cnt=1;mycells=[];
for i=1:ncells
    if ~isempty(datarun.vision.timecourses(i).g) && sum(datarun.vision.timecourses(i).g)~=0
        time_courses(:,cnt)=[datarun.vision.timecourses(i).g(16:end)];        
        isi=diff(datarun.spikes{i}*1000);
        isi=hist(isi,0:1:100);
        isis(:,cnt)=isi(1:end-1);        
        mycells=[mycells cnt];
    end
    cnt=cnt+1;
end


tc_length=length(datarun.vision.timecourses(1).r(16:end));
time_courses=zeros(tc_length*3,ncells);
isis=zeros(100,ncells);
cnt=1;mycells=[];
for i=1:ncells
    if ~isempty(datarun.vision.timecourses(i).r) && sum(datarun.vision.timecourses(i).r)~=0
        time_courses(:,cnt)=[datarun.vision.timecourses(i).r(16:end);...
            datarun.vision.timecourses(i).g(16:end); datarun.vision.timecourses(i).b(16:end)];
        
        isi=diff(datarun.spikes{i}*1000);
        isi=hist(isi,0:1:100);
        isis(:,cnt)=isi(1:end-1);        
        mycells=[mycells cnt];
    end
    cnt=cnt+1;
end



tmp=corr(time_courses);
tmp=tmp-diag(ones(size(tmp,1),1));
close all
mycl={}; cnt=1; thresh=0.98; myMax=1; fl=1; clustered=[];
while fl && cnt<5
    tmp1=tmp;
    tmp1(tmp1<thresh)=0;
    if max(tmp1(:))>0.5
        a=nanmean(tmp1);
        startingCell=find(a==max(a),1);
        
        mycl{cnt}=sort([startingCell find(tmp(startingCell,:)>thresh) ]);
        tmp(mycl{cnt},:)=0;
        tmp(:,mycl{cnt})=0;
        figure(1)
        subplot(4,4,cnt)
        plot(time_courses(:,mycl{cnt}))
        title(int2str(length(mycl{cnt})))
        subplot(4,4,cnt+4)
        plot(isis(:,mycl{cnt}))
        
        subplot(4,4,cnt+8)
        plot_rf_summaries(datarun, datarun.cell_ids(mycl{cnt}), 'clear', false,  'plot_fits', true, 'fit_color', 'k')
        
        clustered=[clustered mycl{cnt}];
        cnt=cnt+1;
    else
        fl=0;
        subplot(4,4,cnt)
        notcl=1:ncells;
        notcl(clustered)=[];
        plot(time_courses(:,notcl))
        title(int2str(length(notcl)))
        
        subplot(4,4,cnt+8)
        plot_rf_summaries(datarun, datarun.cell_ids(mycl{cnt}), 'clear', false,  'plot_fits', true, 'fit_color', 'k')
        
        subplot(4,4,cnt+4)
        plot(isis(:,notcl))
        cnt=cnt+1;
    end
    
end


tmp=corr(isis);
tmp=tmp-diag(ones(size(tmp,1),1));
figure
mycl={}; cnt=1; thresh=0.8; myMax=1; fl=1;
while fl
    tmp1=tmp;
    tmp1(tmp1<thresh)=0;
    if max(tmp1(:))>0.5
        a=nanmean(tmp1);
        startingCell=find(a==max(a),1);
        
        mycl{cnt}=sort([startingCell find(tmp(startingCell,:)>thresh) ]);
        tmp(mycl{cnt},:)=0;
        tmp(:,mycl{cnt})=0;
        subplot(4,4,cnt)
        plot(isis(:,mycl{cnt}))
        title(int2str(length(mycl{cnt})))
        cnt=cnt+1;
    else
        fl=0;
    end
end