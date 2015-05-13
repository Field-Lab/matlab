
% paths

mainpath = '/Volumes/Analysis/2015-03-09-2/';

subpath{1} = 'd02-11-norefit/';
subpath{2} = 'd10-19-norefit/';
subpath{3} = 'd18-28-norefit/';
subpath{4} = 'd27-33-norefit/';
subpath{5} = 'd32-37-norefit/';
subpath{6} = 'd37-40-norefit/';



datarun = load_data([mainpath, subpath{1}, '/data011-from-d02-d11/data011-from-d02-d11']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_ei(datarun,'all');

datarun2 = load_data([mainpath, subpath{2}, '/data011-from-d10-d19/data011-from-d10-d19']);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta','all','keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_neurons(datarun2);
datarun2 = load_ei(datarun2,'all');

datarun3 = load_data([mainpath, subpath{2}, '/data019-from-d10-d19/data019-from-d10-d19']);
datarun3 = load_params(datarun3,'verbose',1);
datarun3 = load_sta(datarun3,'load_sta','all','keep_java_sta',true);
datarun3 = set_polarities(datarun3);
datarun3 = load_neurons(datarun3);
datarun3 = load_ei(datarun3,'all');

datarun4 = load_data([mainpath, subpath{3}, '/data019-from-d18-28/data019-from-d18-28']);
datarun4 = load_params(datarun4,'verbose',1);
datarun4 = load_sta(datarun4,'load_sta','all','keep_java_sta',true);
datarun4 = set_polarities(datarun4);
datarun4 = load_neurons(datarun4);
datarun4 = load_ei(datarun4,'all');


[a,b]=map_ei(datarun,datarun2);

[c,d]=map_ei(datarun,datarun4);


filepath='/Users/alexth/Desktop/Light_adaptation//2015-03-09-2/cell_mapping/';

for i=1:length(a)
    
    
    for j=1:length(datarun.cell_types)
        if ~isempty(find(datarun.cell_ids(i)==datarun.cell_types{j}.cell_ids, 1))
            folder=datarun.cell_types{j}.name;
            if ~exist([filepath,folder],'dir')
                mkdir([filepath,folder]);
            end
            break
        end
    end
    
    fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[201         239        1183         860]);
    
    sta = squeeze(datarun.stas.stas{i});
    subplot(2,2,1)
    colormap gray
    [t,ind]=max(abs(sta(:)));
    [~,~,t]=ind2sub(size(sta),ind);
    imagesc(sta(:,:,t));
    set(gca,'dataaspectratio',[1 1 1])
    title(['Cell ', int2str(datarun.cell_ids(i)), ', data011 from 02-11'])
    
    if ~isempty(a{i})
        cellID2=find(datarun2.cell_ids==a{i});
        sta = squeeze(datarun2.stas.stas{cellID2});
        subplot(2,2,2)
        colormap gray
        [t,ind]=max(abs(sta(:)));
        [~,~,t]=ind2sub(size(sta),ind);
        imagesc(sta(:,:,t));
        set(gca,'dataaspectratio',[1 1 1])
        title(['Cell ', int2str(a{i}), ', data011 from 10-19'])
        
        sta = squeeze(datarun3.stas.stas{cellID2});
        subplot(2,2,3)
        colormap gray
        [t,ind]=max(abs(sta(:)));
        [~,~,~,p]=ind2sub(size(sta),ind);
        sta=sta/(2*t)+0.5;
        imagesc(sta(:,:,:,p));
        set(gca,'dataaspectratio',[1 1 1])
        title(['Cell ', int2str(a{i}), ', data019 from 10-19'])
    end
    
    
    if ~isempty(c{i})
        cellID3=find(datarun4.cell_ids==c{i});
        sta = squeeze(datarun4.stas.stas{cellID3});
        subplot(2,2,4)
        colormap gray
        [t,ind]=max(abs(sta(:)));
        [~,~,~,p]=ind2sub(size(sta),ind);
        sta=sta/(2*t)+0.5;
        imagesc(sta(:,:,:,p));
        set(gca,'dataaspectratio',[1 1 1])
        title(['Cell ', int2str(c{i}), ', data019 from 18-28'])
    end
    
    
    
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(datarun.cell_ids(i))]));
    close(fig)
end
    
    




datarun = load_data([mainpath, subpath{6}, '/data037-from-d37-40/data037-from-d37-40']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_ei(datarun,'all');

datarun2 = load_data([mainpath, subpath{2}, '/data015-from-d10-d19/data015-from-d10-d19']);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta','all','keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_neurons(datarun2);
datarun2 = load_ei(datarun2,'all');

[a,b]=map_ei(datarun,datarun2);


filepath='/Users/alexth/Desktop/Light_adaptation//2015-03-09-2/cell_mapping/d15-d37/';

for i=1:length(a)
    
    
    for j=1:length(datarun.cell_types)
        if ~isempty(find(datarun.cell_ids(i)==datarun.cell_types{j}.cell_ids, 1))
            folder=datarun.cell_types{j}.name;
            if ~exist([filepath,folder],'dir')
                mkdir([filepath,folder]);
            end
            break
        end
    end
    
    fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[201         239        1183         860]);
    
    sta = squeeze(datarun.stas.stas{i});
    subplot(1,2,1)
    colormap gray
    [t,ind]=max(abs(sta(:)));
    [~,~,t]=ind2sub(size(sta),ind);
    imagesc(sta(:,:,t));
    set(gca,'dataaspectratio',[1 1 1])
    title(['Cell ', int2str(datarun.cell_ids(i)), ', data037 from 37-40'])
    
    if ~isempty(a{i})
        cellID2=find(datarun2.cell_ids==a{i});
        sta = squeeze(datarun2.stas.stas{cellID2});
        subplot(1,2,2)
        colormap gray
        [t,ind]=max(abs(sta(:)));
        [~,~,t]=ind2sub(size(sta),ind);
        imagesc(sta(:,:,t));
        set(gca,'dataaspectratio',[1 1 1])
        title(['Cell ', int2str(a{i}), ', data015 from 10-19'])
        
    end
    
    
    
    
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(datarun.cell_ids(i))]));
    close(fig)
end
    
    


datarun = load_data([mainpath, subpath{6}, '/data038-from-d37-40/data038-from-d37-40']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_ei(datarun,'all');

datarun2 = load_data([mainpath, subpath{1}, '/data002-from-d02-d11/data002-from-d02-d11']);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2,'load_sta','all','keep_java_sta',true);
datarun2 = set_polarities(datarun2);
datarun2 = load_neurons(datarun2);
datarun2 = load_ei(datarun2,'all');

figure
plot_rf_summaries(datarun, {5}, 'clear', false, 'scale',2, 'plot_fits', true, 'fit_color', 'r')
plot_rf_summaries(datarun2, {5}, 'clear', false,  'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun, {3}, 'clear', false,  'plot_fits', true, 'fit_color', [0 0.6 0])
plot_rf_summaries(datarun, {4}, 'clear', false,  'plot_fits', true, 'fit_color', [0 0.8 0.8])

