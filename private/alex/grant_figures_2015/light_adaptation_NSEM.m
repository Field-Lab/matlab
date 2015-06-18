clear


load('/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/data_011_015_chopped.mat','data1','data2')

data='015';
starun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',data,'-from-d05-d27/data',data,'-from-d05-d27']);
starun = load_params(starun,'verbose',1);
filepath = '/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/nd3_2_chopped/';
nds={'NDF3', 'NDF2'};


sr=120;
sig=50;
st=10000/sr*6.2/(60*sig);
time_list=-3.1:st:3.1;
kern=zeros(1,length(time_list));
for i=1:length(time_list)    
    kern(i)=250/sig*exp((1-time_list(i)^2)/2);
end

figure
plot(kern)

0.8/8.24


%% OFF midget, cell 2945
all_time = 1:1/120:30;
my_times = [1, 12,19, 26]; % in s
for j=1:length(my_times)
    mt(j) = length(1:1/120:my_times(j));
end
for cellID = 2945
    
    datarunID=find(starun.cell_ids==cellID);
    
    %%%%%% NDF2 raster
    figsize = [260 40];
    figure
    set(gcf,'units', 'points','position',[400   400  figsize(1), figsize(2)])    
    h = subplot('position',[0 0 1 1]);
    
    spikes=data2.nm.raw_spikes{datarunID};
    trigs = data2.nm.trigs;
    myTrigs=[0 find(diff(trigs)>0.9)'];
    splitSpikes=cell(14,1);
    for j=6:19
        tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
            - trigs(myTrigs(j)+1);
        tmp = [ tmp(tmp>1&tmp<12); tmp(tmp>19&tmp<26)-(19-12)+1];
        splitSpikes{j-5}=tmp*1000;
    end    
    
    tt = 0;
    for j=1:14
        tt = tt +sum(splitSpikes{j}>7000 & splitSpikes{j}<8000);
    end
    tt = tt/14;

    mean(tmp(length(1:1/120:7):length(1:1/120:8)))
    
    hold on
    for j=1:14
        plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
    end
    
%     for j = 1:5
%         line([12, 13]*1000, [0, 0.2]+0.2*j-0.2, 'linestyle', '-','color', [1 1 1]*0.3)
%     end
    
    axis([1000 (26-(19-12)+1)*1000, -0.02 0.8])
    set(gca, 'visible', 'off')

    %%%%%% NDF2 FR

    figsize = [260 50];
    figure
    set(gcf,'units', 'points','position',[400   400  figsize(1), figsize(2)])    
    subplot('position',[0 0 1 1])
    hold on    
  
    % actual
    tt = unique(data2.nm.asr{datarunID});
    tmp = data2.nm.asr{datarunID}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(tmp(mt(1):mt(2)),'k', 'linewidth',2) 
    plot((mt(2)+120):(mt(2)+mt(4)-mt(3)+120), tmp(mt(3):mt(4)),'k', 'linewidth',2) 
    
    % predicted
    tt = unique(data2.nm.psr{datarunID});
    tmp = data2.nm.psr{datarunID}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(tmp(mt(1):mt(2)),'r', 'linewidth',1) 
    hold on
    plot((mt(2)+120):(mt(2)+mt(4)-mt(3)+120), tmp(mt(3):mt(4)),'r', 'linewidth',1)  

%     for j = 1:5
%         line([mt(2), mt(2)+120], [0, 0.2]+0.2*j-0.2, 'linestyle', '-','color', [1 1 1]*0.3)
%     end
    
%     set(gca,'xtick', [mt(1:2), (mt(2)+100),(mt(2)+mt(4)-mt(3)+100) ], 'xticklabel', {'1', '12', '19', '26'})
    axis([0 2262, 0 0.8])
    set(gca, 'visible', 'off')
    
    %%%%%% NDF3 FR
    figsize = [260 50];
    figure
    set(gcf,'units', 'points','position',[400   400  figsize(1), figsize(2)])    
    subplot('position',[0 0 1 1])
    hold on    
  
    % actual
    tt = unique(data1.nm.asr{datarunID});
    tmp = data1.nm.asr{datarunID}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(tmp(mt(1):mt(2)),'k', 'linewidth',2) 
    plot((mt(2)+120):(mt(2)+mt(4)-mt(3)+120), tmp(mt(3):mt(4)),'k', 'linewidth',2) 
    
    % predicted
    tt = unique(data1.nm.psr{datarunID});
    tmp = data1.nm.psr{datarunID}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(tmp(mt(1):mt(2)),'r', 'linewidth',1) 
    hold on
    plot((mt(2)+120):(mt(2)+mt(4)-mt(3)+120), tmp(mt(3):mt(4)),'r', 'linewidth',1)  

%     for j = 1:5
%         line([mt(2), mt(2)+120], [0, 0.2]+0.2*j-0.2, 'linestyle', '-','color', [1 1 1]*0.3)
%     end
%     
%     set(gca,'xtick', [mt(1:2), (mt(2)+100),(mt(2)+mt(4)-mt(3)+100) ], 'xticklabel', {'1', '12', '19', '26'})
    axis([0 2262, 0 0.8])
    set(gca, 'visible', 'off')

    
    %%%%%% NDF3 raster
    
    figsize = [260 40];
    figure
    set(gcf,'units', 'points','position',[400   400  figsize(1), figsize(2)])    
    h = subplot('position',[0 0 1 1]);
    
    spikes=data1.nm.raw_spikes{datarunID};
    trigs = data1.nm.trigs;
    myTrigs=[0 find(diff(trigs)>0.9)'];
    splitSpikes=cell(14,1);
    for j=6:19
        tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
            - trigs(myTrigs(j)+1);
        tmp = [ tmp(tmp>1&tmp<12); tmp(tmp>19&tmp<26)-(19-12)+1];
        splitSpikes{j-5}=tmp*1000;
    end    
    
    hold on
    for j=1:14
        plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
    end
    
%     for j = 1:5
%         line([12, 13]*1000, [0, 0.2]+0.2*j-0.2, 'linestyle', '-','color', [1 1 1]*0.3)
%     end
    
    axis([1000 (26-(19-12)+1)*1000, -0.02 0.8])
    set(gca, 'visible', 'off')
    
    
end



 %% ON parasol, cell 2089

all_time = 1:1/120:30;
my_times = [1, 22]; % in s
clear mt
for j=1:length(my_times)
    mt(j) = length(1:1/120:my_times(j));
end

for cellID = 2089
    
    datarunID=find(starun.cell_ids==cellID);
    
    %%%%%% NDF2 raster
    figsize = [260 40];
    figure
    set(gcf,'units', 'points','position',[400   400  figsize(1), figsize(2)])    
    h = subplot('position',[0 0 1 1]);
    
    spikes=data2.nm.raw_spikes{datarunID};
    trigs = data2.nm.trigs;
    myTrigs=[0 find(diff(trigs)>0.9)'];
    splitSpikes=cell(14,1);
    for j=6:19
        tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
            - trigs(myTrigs(j)+1);
        tmp = [ tmp(tmp>1&tmp<22)];
        splitSpikes{j-5}=tmp*1000;
    end    
    
    
        
    tt = 0;
    for j=1:14
        tt = tt +sum(splitSpikes{j}>8000 & splitSpikes{j}<9000);
    end
    tt = tt/14;

    mean(tmp(length(1:1/120:8):length(1:1/120:9)))
    
    hold on
    for j=1:14
        plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+3/13*j-3/13, '.k', 'markersize',0.1)
    end
    
    axis([1000 22*1000, -0.02 3])
    set(gca, 'visible', 'off')


    %%%%%% NDF2 FR

    
    figsize = [260 50];
    figure
    set(gcf,'units', 'points','position',[400   400  figsize(1), figsize(2)])    
    subplot('position',[0 0 1 1])
    hold on    
  
    % actual
    tt = unique(data2.nm.asr{datarunID});
    tmp = data2.nm.asr{datarunID}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(tmp(mt(1):mt(2)),'k', 'linewidth',2) 
    
    % predicted
    tt = unique(data2.nm.psr{datarunID});
    tmp = data2.nm.psr{datarunID}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(tmp(mt(1):mt(2)),'r', 'linewidth',1)  

    axis([0 mt(2)+1, 0 3])
    set(gca, 'visible', 'off')
        

    
    %%%%%% NDF3 FR
    figsize = [260 50];
    figure
    set(gcf,'units', 'points','position',[400   400  figsize(1), figsize(2)])    
    subplot('position',[0 0 1 1])
    hold on    
  
    % actual
    tt = unique(data1.nm.asr{datarunID});
    tmp = data1.nm.asr{datarunID}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(tmp(mt(1):mt(2)),'k', 'linewidth',2) 
    
    % predicted
    tt = unique(data1.nm.psr{datarunID});
    tmp = data1.nm.psr{datarunID}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(tmp(mt(1):mt(2)),'r', 'linewidth',1)  

    axis([0 mt(2)+1, 0 3])
    set(gca, 'visible', 'off')
    

    
    %%%%%% NDF3 raster
    %%%%%% NDF2 raster
    figsize = [260 40];
    figure
    set(gcf,'units', 'points','position',[400   400  figsize(1), figsize(2)])    
    h = subplot('position',[0 0 1 1]);
    
    spikes=data1.nm.raw_spikes{datarunID};
    trigs = data1.nm.trigs;
    myTrigs=[0 find(diff(trigs)>0.9)'];
    splitSpikes=cell(14,1);
    for j=6:19
        tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
            - trigs(myTrigs(j)+1);
        tmp = [ tmp(tmp>1&tmp<22)];
        splitSpikes{j-5}=tmp*1000;
    end    
    
    hold on
    for j=1:14
        plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+3/13*j-3/13, '.k', 'markersize',0.1)
    end
    
    axis([1000 22*1000, -0.02 3])
    set(gca, 'visible', 'off')
    
    
end
