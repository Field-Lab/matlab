
clear

%%%% 2012-09-13-2 INFO
onM_d01 = [1909 2225 2267 2449 2795 4056 4337 5556 6079 6153 6409 6679 6802 6917 7292];

offM_d05 = [1396 1607 1727 1816 1966 2026 2116 2386 2942 6991];
offP_d05 = [46 376 1231 4441 5371 5551 6496];

offM_d09 = [151 1396 1741 1816 1863 1966 2026 2116 6796 6991];
offP_d09 = [376 1471 4636 4771 5551];


piece = '2012-09-13-2';
run = 'data013';

%%%% 2013-10-10-5 INFO %%%%%%

offM{1}=[721 726 2809 3005 3183 4222 4758 5090] % from data002 for data003 for 2013-10-10-5, 117 cones, stimulated only 92?

piece = '2013-10-10-5';
run = 'data003';

% define data path

datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);

datarun=load_data(datarun,opt);
    
filepath='/Volumes/Analysis/';


% stimulus=read_stim_lisp_output_ath('2012-09-13-2','s13','allconesd09/')
stimulus=read_stim_lisp_output_ath('2013-10-10-5','s03')

parsed=parse_stim_rgbs_ath(stimulus);

load('/Volumes/Analysis/2012-09-13-2/stimuli/allconesd09/map-0000.txt')
figure
imagesc(map_0000)
length(unique(map_0000))

load(['/Volumes/Acquisition/maps/2013-10-10-5/allconesOffmidgets/map-0001.txt'])
figure
imagesc(map_0001)
length(unique(map_0001))


a=parsed.rgbs;

k=cell(1,stimulus.numcones); % stimuli list for each cone
m=zeros(size(stimulus.pulses,2)/stimulus.numcones,stimulus.numcones); % when stimuli were applied
for i=1:size(stimulus.pulses,2)
    p=a{i};
    t=find(p(:,1));
    k{t}=[k{t} p(t,1)];
    p=find(m(:,t)==0,1);
    m(p,t)=i;
end





tr=datarun.triggers(1:2:end);
mmin=0;mmax=0.75;

for i=10:15%408
    spikes=datarun.spikes{i};
    conv=zeros(29,750);
    figure
    for j=1:29
        trig_local=m(k{j}==-0.48,j);
        for f=trig_local'
            h=spikes-tr(f);
            hh= h>=mmin & h<mmax;
            hh=round(h(hh)*1000);
            tmp=convolved(hh,40,750);
            tmp=tmp(121:end-120);
            conv(j,:)=conv(j,:)+tmp;
        end
        subplot(5,6,j)
        plot(conv(j,:)/50)
        hold on
        
        conv(j,:)=conv(j,:)-conv(j,:);
        trig_local=m(k{j}==-0.288,j);
        for f=trig_local'
            h=spikes-tr(f);
            hh= h>=mmin & h<mmax;
            hh=round(h(hh)*1000);
            tmp=convolved(hh,40,750);
            tmp=tmp(121:end-120);
            conv(j,:)=conv(j,:)+tmp;
        end
        plot(conv(j,:)/50,'r')
        
        axis([0 750 0 30])
    end
  
end




piece = '2012-09-13-2';
run = 'data009';

% define data path
datarunA = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);

datarunA = load_params(datarunA,struct('verbose',1));

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);

datarunA=load_data(datarunA,opt);


[cell_list, failed_cells]=map_ei(datarun,datarunA);
cellOfInterest=4486;
a=0;
for i=1:408
    if ~isempty(cell_list{i})
        a=a+1;
    end
    if isequal(cell_list{i},cellOfInterest)
        indInDatarun=i;
    end
end
cellIndInDatarun=datarun.cell_ids(indInDatarun);

kmaps=map_0000;
for i=indInDatarun%408
    spikes=datarun.spikes{i};
    conv=zeros(29,750);
    figure
    for j=1:29
        trig_local=m(k{j}==-0.48,j);
        for f=trig_local'
            h=spikes-tr(f);
            hh= h>=mmin & h<mmax;
            hh=round(h(hh)*1000);
            tmp=convolved(hh,40,750);
            tmp=tmp(121:end-120);
            conv(j,:)=conv(j,:)+tmp;
        end
        subplot(5,6,j)
        plot(conv(j,:)/50)
        kkk=conv(j,:)/50;
        hold on
        
        conv(j,:)=conv(j,:)-conv(j,:);
        trig_local=m(k{j}==-0.288,j);
        for f=trig_local'
            h=spikes-tr(f);
            hh= h>=mmin & h<mmax;
            hh=round(h(hh)*1000);
            tmp=convolved(hh,40,750);
            tmp=tmp(121:end-120);
            conv(j,:)=conv(j,:)+tmp;
        end
        plot(conv(j,:)/50,'r')
        
        kmaps(map_0000==j)=max(kkk(1:300))-mean(kkk(450:650));
    end
  
end

figure
imagesc(kmaps)


for i=indInDatarun%408
    figure
    for j=1:29
        subplot(5,6,j)
        kmaps=map_0000;
        kmaps(map_0000~=j)=0;
        imagesc(kmaps)
    end
  
end

for i=indInDatarun%408
    for j=1:29
%         subplot(5,6,j)                                    
figure
        kmaps=map_0000;
        kmaps(map_0000~=j)=0;
        imagesc(kmaps)
    end
  
end

%%%%%%%% datarunA cone plotting %%%%%%%%
% load data
datarunA = load_params(datarunA,struct('verbose',1));  
datarunA = load_sta(datarunA,'load_sta',[]);
datarunA = set_polarities(datarunA);
% if error loading cones - probably multiple cone maps saved, put the number as second argument
datarunA = load_cones(datarunA, 'standard'); % load_cones(datarun, 'Analysis');
datarunA = make_mosaic_struct(datarunA);
datarunA = get_sta_fits_from_vision(datarunA);  
figure

figure
imagesc(map_0000)
hold on
plot_rf_summaries(datarunA, {4}, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
plot_rf_summaries(datarunA, {2}, 'scale',2,'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')

offM_d09 = [151 1396 1741 1816 1863 1966 2026 2116 6796 6991];
offP_d09 = [376 1471 4636 4771 5551];

for i=1:size(datarunA.cell_types,2)
   a=datarunA.cell_types{1, i}.cell_ids;
   intersect(a,offP_d09)
end



figure
plot_rf_summaries(datarunA, {2}, 'clear', false, 'label', true, 'label_color', 'k', 'plot_fits', true, 'fit_color', 'k')
hold on
plot_rf_summaries(datarun, {2}, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')


figure
imagesc(map_0000)
hold on
plot_rf_summaries(datarunA, {2}, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
plot_rf_summaries(datarun, {2}, 'scale',2,'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')