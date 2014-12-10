%% Setup
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
staopts = loadopts;
staopts.load_sta = true;
staopts.load_sta_params = struct('load_sta', [], 'guess_stimulus', false);
staopts.set_polarities = {'guess', true};
pieces = containers.Map('KeyType', 'char', 'ValueType', 'any');
keep_vars = {'pieces'; 'loadopts'; 'staopts'; 'keep_vars'};

%% Load 2011-12-13-2
piece = '2011-12-13-2';

% d08s = load_data([piece '/streamed/data008-0/data008-0'], staopts);
% d08s = load_cones_ath(d08s, 'bayes');
% d08s = conepreprocess_wnm(d08s);
% conepreprocess_save(d08s,'date', piece);

load(['/Volumes/Analysis/', piece, '/subunits/data008-0/conepreprocess.mat']);
d08s = datarun;
clear datarun myrgc

rgc = [1321 1351 2251 3586 4576 5162];
for i=1:length(rgc)
    myrgc(i)=find(d08s.cell_ids==rgc(i));
end

d08s.channels=d08s.channels(myrgc);
d08s.cell_ids=d08s.cell_ids(myrgc);
d08s.ei.eis=d08s.ei.eis(myrgc);
d08s.ei.num_spikes=d08s.ei.num_spikes(myrgc);
d08s.stas.rfs=d08s.stas.rfs(myrgc);
d08s.stas.marks=d08s.stas.marks(myrgc);
d08s.stas.rf_coms=d08s.stas.rf_coms(myrgc);
d08s.stas.time_courses=d08s.stas.time_courses(myrgc);
d08s.stas.fits=d08s.stas.fits(myrgc);
d08s.stas.polarities=d08s.stas.polarities(myrgc);
d08s.cones.weights=d08s.cones.weights(:, myrgc);
d08s.spike_rate=d08s.spike_rate(myrgc, :);

d10 = load_data([piece '/data010'], loadopts);
d10 = read_stim_lisp_output(d10, ':2011-12-13-2_f08_allcones');
d10.stimulus = convert_stimulus_to_combined_maps(d10.stimulus);
d10.stimulus = parse_stim_rgbs(d10.stimulus);
d10.stimulus.triggers = d10.triggers(1:2:end);
d10.mapd08s = map_ei(d08s, d10);
raster_rgcs=[];
for rgcs=rgc    
    tmp = d10.mapd08s{get_cell_indices(d08s, rgcs)};
    raster_rgcs=[raster_rgcs find(d10.cell_ids==tmp)];
end

d10.channels=d10.channels(raster_rgcs);
d10.spikes=d10.spikes(raster_rgcs);
d10.cell_ids=d10.cell_ids(raster_rgcs);
d10.ei.eis=d10.ei.eis(raster_rgcs);
d10.ei.num_spikes=d10.ei.num_spikes(raster_rgcs);

savepath='/Volumes/Analysis/2011-12-13-2/singlecones/';
mkdir(savepath)
save([savepath, 'dataruns.mat'], 'd08s','d10');



%% Load 2012-04-13-1
piece = '2012-04-13-1';

% d02s = load_data([piece '/streamed/data002/data002'], staopts);
% d02s = load_cones_ath(d02s, 'acquisition');
% d02s.names.rrs_movie_path='/Volumes/Analysis/2012-04-13-1/streamed/data002/data002.movie';
% d02s = conepreprocess_wnm(d02s);
% conepreprocess_save(d02s,'date', piece);

% d06s = load_data([piece '/streamed/data006/data006'], staopts);
% d06s = load_cones_ath(d06s, piece);
% d06s.names.rrs_movie_path='/Volumes/Analysis/2012-04-13-1/streamed/data006/data006.movie';
% d06s = conepreprocess_wnm(d06s);
% conepreprocess_save(d06s,'date', piece);


load(['/Volumes/Analysis/', piece, '/subunits/data002/conepreprocess.mat']);
d02s = datarun;
clear datarun myrgc

rgc = 2389;
for i=1:length(rgc)
    myrgc(i)=find(d02s.cell_ids==rgc(i));
end

d02s.channels=d02s.channels(myrgc);
d02s.cell_ids=d02s.cell_ids(myrgc);
d02s.ei.eis=d02s.ei.eis(myrgc);
d02s.ei.num_spikes=d02s.ei.num_spikes(myrgc);
d02s.stas.rfs=d02s.stas.rfs(myrgc);
d02s.stas.marks=d02s.stas.marks(myrgc);
d02s.stas.rf_coms=d02s.stas.rf_coms(myrgc);
d02s.stas.time_courses=d02s.stas.time_courses(myrgc);
d02s.stas.fits=d02s.stas.fits(myrgc);
d02s.stas.polarities=d02s.stas.polarities(myrgc);
d02s.cones.weights=d02s.cones.weights(:, myrgc);
d02s.spike_rate=d02s.spike_rate(myrgc, :);


load(['/Volumes/Analysis/', piece, '/subunits/data006/conepreprocess.mat']);
d06s = datarun;
clear datarun myrgc

rgc = 6136;
for i=1:length(rgc)
    myrgc(i)=find(d06s.cell_ids==rgc(i));
end

d06s.channels=d06s.channels(myrgc);
d06s.cell_ids=d06s.cell_ids(myrgc);
d06s.ei.eis=d06s.ei.eis(myrgc);
d06s.ei.num_spikes=d06s.ei.num_spikes(myrgc);
d06s.stas.rfs=d06s.stas.rfs(myrgc);
d06s.stas.marks=d06s.stas.marks(myrgc);
d06s.stas.rf_coms=d06s.stas.rf_coms(myrgc);
d06s.stas.time_courses=d06s.stas.time_courses(myrgc);
d06s.stas.fits=d06s.stas.fits(myrgc);
d06s.stas.polarities=d06s.stas.polarities(myrgc);
d06s.cones.weights=d06s.cones.weights(:, myrgc);
d06s.spike_rate=d06s.spike_rate(myrgc, :);


d09cm = load_data([piece '/d00-01-02-04-05-06-07-09-10-norefit/data009-from-data000_data001_data002_data004_data005_data006_data007_data009_data010/data009-from-data000_data001_data002_data004_data005_data006_data007_data009_data010'], loadopts);
d09cm = read_stim_lisp_output(d09cm, ':2012-04-13-1:f06_allcones');
d09cm.stimulus = convert_stimulus_to_combined_maps(d09cm.stimulus);
d09cm.stimulus = parse_stim_rgbs(d09cm.stimulus);
d09cm.stimulus.triggers = d09cm.triggers(1:2:end);
d09cm.mapd06s = map_ei(d06s, d09cm);

raster_rgcs=[];
for rgcs=rgc    
    tmp = d09cm.mapd06s{get_cell_indices(d06s, rgcs)};
    raster_rgcs=[raster_rgcs find(d09cm.cell_ids==tmp)];
end

d09cm.channels=d09cm.channels(raster_rgcs);
d09cm.spikes=d09cm.spikes(raster_rgcs);
d09cm.cell_ids=d09cm.cell_ids(raster_rgcs);
d09cm.ei.eis=d09cm.ei.eis(raster_rgcs);
d09cm.ei.num_spikes=d09cm.ei.num_spikes(raster_rgcs);

savepath='/Volumes/Analysis/2012-04-13-1/singlecones/';
mkdir(savepath)
save([savepath, 'dataruns.mat'], 'd02s','d06s','d09cm');



%% Load 2012-09-06-0
piece = '2012-09-06-0';
% 
% d04s = load_data([piece, '/streamed/data004/data004'], staopts);
% d04s = load_cones_ath(d04s, piece);
% d04s.names.rrs_movie_path='/Volumes/Analysis/2012-09-06-0/streamed/data004/data004.movie';
% d04s = conepreprocess_wnm(d04s);
% conepreprocess_save(d04s,'date', piece, 'cell_types', {10});

load(['/Volumes/Analysis/', piece, '/subunits/data004/conepreprocess.mat']);
d04s = datarun;
clear datarun myrgc

rgc = [5103 5746 6031 6317 7022];
for i=1:length(rgc)
    myrgc(i)=find(d04s.cell_ids==rgc(i));
end

d04s.channels=d04s.channels(myrgc);
d04s.cell_ids=d04s.cell_ids(myrgc);
d04s.ei.eis=d04s.ei.eis(myrgc);
d04s.ei.num_spikes=d04s.ei.num_spikes(myrgc);
d04s.stas.rfs=d04s.stas.rfs(myrgc);
d04s.stas.marks=d04s.stas.marks(myrgc);
d04s.stas.rf_coms=d04s.stas.rf_coms(myrgc);
d04s.stas.time_courses=d04s.stas.time_courses(myrgc);
d04s.stas.fits=d04s.stas.fits(myrgc);
d04s.stas.polarities=d04s.stas.polarities(myrgc);
d04s.cones.weights=d04s.cones.weights(:, myrgc);
d04s.spike_rate=d04s.spike_rate(myrgc, :);


d08 = load_data([piece '/data008'], loadopts);
d08 = read_stim_lisp_output(d08);
d08.stimulus = parse_stim_rgbs(d08.stimulus);
d08.stimulus.triggers = d08.triggers(1:2:end);
d08.mapd04s = map_ei(d04s, d08);

raster_rgcs=[];
for rgcs=rgc    
    tmp = d08.mapd04s{get_cell_indices(d04s, rgcs)};
    raster_rgcs=[raster_rgcs find(d08.cell_ids==tmp)];
end

d08.channels=d08.channels(raster_rgcs);
d08.spikes=d08.spikes(raster_rgcs);
d08.cell_ids=d08.cell_ids(raster_rgcs);
d08.ei.eis=d08.ei.eis(raster_rgcs);
d08.ei.num_spikes=d08.ei.num_spikes(raster_rgcs);

savepath='/Volumes/Analysis/2012-09-06-0/singlecones/';
mkdir(savepath)
save([savepath, 'dataruns.mat'], 'd04s','d08');


%% Load 2012-09-21-2
piece = '2012-09-21-2';

% d09 = load_data([piece '/data009'], staopts);
% d09 = load_cones_ath(d09, piece);
% d09 = conepreprocess_wnm(d09);
% conepreprocess_save(d09,'date', piece, 'cell_types', {8});

load(['/Volumes/Analysis/', piece, '/subunits/data009/conepreprocess.mat']);
d09 = datarun;
clear datarun myrgc

rgc = [5086  6346 6406 7696];
for i=1:length(rgc)
    myrgc(i)=find(d09.cell_ids==rgc(i));
end

d09.channels=d09.channels(myrgc);
d09.cell_ids=d09.cell_ids(myrgc);
d09.ei.eis=d09.ei.eis(myrgc);
d09.ei.num_spikes=d09.ei.num_spikes(myrgc);
d09.stas.rfs=d09.stas.rfs(myrgc);
d09.stas.marks=d09.stas.marks(myrgc);
d09.stas.rf_coms=d09.stas.rf_coms(myrgc);
d09.stas.time_courses=d09.stas.time_courses(myrgc);
d09.stas.fits=d09.stas.fits(myrgc);
d09.stas.polarities=d09.stas.polarities(myrgc);
d09.cones.weights=d09.cones.weights(:, myrgc);
d09.spike_rate=d09.spike_rate(myrgc, :);


d13 = load_data([piece '/data013'], loadopts);
d13 = read_stim_lisp_output(d13);
d13.stimulus = parse_stim_rgbs(d13.stimulus);
d13.stimulus.triggers = d13.triggers(1:2:end);
d13.mapd09 = map_ei(d09, d13);


raster_rgcs=[];
for rgcs=rgc    
    tmp = d13.mapd09{get_cell_indices(d09, rgcs)};
    raster_rgcs=[raster_rgcs find(d13.cell_ids==tmp)];
end

d13.channels=d13.channels(raster_rgcs);
d13.spikes=d13.spikes(raster_rgcs);
d13.cell_ids=d13.cell_ids(raster_rgcs);
d13.ei.eis=d13.ei.eis(raster_rgcs);
d13.ei.num_spikes=d13.ei.num_spikes(raster_rgcs);

savepath='/Volumes/Analysis/2012-09-21-2/singlecones/';
mkdir(savepath)
save([savepath, 'dataruns.mat'], 'd09','d13');

%% Load 2012-08-21-2
piece = '2012-08-21-2';

% d01s = load_data([piece '/streamed/data001/data001'], staopts);
% d01s = load_cones_ath(d01s, 2);
% d01s = conepreprocess_wnm(d01s);
% conepreprocess_save(d01s,'date', piece);

load(['/Volumes/Analysis/', piece, '/subunits/data001/conepreprocess.mat']);
d01s = datarun;
clear datarun myrgc

rgc = [391 2266 2828 5447 7682];
for i=1:length(rgc)
    myrgc(i)=find(d01s.cell_ids==rgc(i));
end

d01s.channels=d01s.channels(myrgc);
d01s.cell_ids=d01s.cell_ids(myrgc);
d01s.ei.eis=d01s.ei.eis(myrgc);
d01s.ei.num_spikes=d01s.ei.num_spikes(myrgc);
d01s.stas.rfs=d01s.stas.rfs(myrgc);
d01s.stas.marks=d01s.stas.marks(myrgc);
d01s.stas.rf_coms=d01s.stas.rf_coms(myrgc);
d01s.stas.time_courses=d01s.stas.time_courses(myrgc);
d01s.stas.fits=d01s.stas.fits(myrgc);
d01s.stas.polarities=d01s.stas.polarities(myrgc);
d01s.cones.weights=d01s.cones.weights(:, myrgc);
d01s.spike_rate=d01s.spike_rate(myrgc, :);
d01s.cones.rf_fits=d01s.cones.rf_fits(myrgc);


d03 = load_data([piece '/data003'], loadopts);
d03 = read_stim_lisp_output(d03);
d03.stimulus = parse_stim_rgbs(d03.stimulus);
d03.stimulus.triggers = d03.triggers(1:2:end);
d03.mapd01s = map_ei(d01s, d03);


raster_rgcs=[];
for rgcs=rgc    
    tmp = d03.mapd01s{get_cell_indices(d01s, rgcs)};
    raster_rgcs=[raster_rgcs find(d03.cell_ids==tmp)];
end

d03.channels=d03.channels(raster_rgcs);
d03.spikes=d03.spikes(raster_rgcs);
d03.cell_ids=d03.cell_ids(raster_rgcs);
d03.ei.eis=d03.ei.eis(raster_rgcs);
d03.ei.num_spikes=d03.ei.num_spikes(raster_rgcs);

savepath='/Volumes/Analysis/2012-08-21-2/singlecones/';
mkdir(savepath)
save([savepath, 'dataruns.mat'], 'd01s','d03');

%% load all data and do Peter's stuff

dates = ['2011-12-13-2'; '2012-04-13-1'; '2012-09-06-0'; '2012-09-21-2'; '2012-08-21-2'];
coneruns = {'d08s', 'd06s', 'd04s', 'd09', 'd01s'};
rasterruns = {'d10', 'd09cm', 'd08', 'd13', 'd03'};
ncells=[1 1 1 1 1 1 2 3 3 3 3 3 4 4 4 4 5 5 5 5 5];
rgcs = [1321 1351 2251 3586 4576 5162 6136 5103 5746 6031 6317 7022 5086 6346 6406 7696 391 2266 2828 5447 7682];
padfactors  = [6 6 6 6 6 6 5 4 4 5 4 7 5 4.2 8 5 5.8 5 4 5 1.8];
badregionss = {[] [10] [8 9] [13] [6 7 8 10 11] [] [] [] [] [] [] [] [] [16] [15] [] [] [] [] [] []};
ars = [1 1 1 1 1 1 1 1 1 1 1 1 1 1.2 1 1 1 1 1 1 1];
for i=1:size(coneruns,2)
    load(['/Volumes/Analysis/', dates(i,:), '/singlecones/dataruns.mat'])
end


sepfigs = false;
if ~sepfigs, figure; end
clear rfweights crweights f gof

coneID=cell(1, length(ncells));

for i = 1:length(ncells)
    
    conerun = eval(coneruns{ncells(i)});
    rasterrun = eval(rasterruns{ncells(i)});
    rasterrun.map=eval(['rasterrun.map' coneruns{ncells(i)}]);
    conerun.rgcs = rgcs(i);
    
    rasterrun.rgcs = rasterrun.map(get_cell_indices(conerun, conerun.rgcs));
    
    padfactor = padfactors(i);
    badregions = badregionss{i};
    ar = ars(i);
   
    localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
    localmap = localmask.*rasterrun.stimulus.mapims{1};
    if ~isempty(badregions)
        for r = badregions, localmap(localmap == r) = 0; end
    end
    
    bb=full(localmap);
    aa=sub2ind(size(localmap),round(conerun.cones.centers(:,1)*2),round(conerun.cones.centers(:,2)*2));
    un=unique(bb)';
    cnt=1;
    for j=un(2:end)
        r=find(bb'==j);
        conetmp=intersect(aa, r);
        [c1,c2]=ind2sub(size(localmap),conetmp);
        coneID{i}(cnt)=find(round(conerun.cones.centers(:,1)*2)==c1 & round(conerun.cones.centers(:,2)*2)==c2);
        cnt=cnt+1;
    end
    
    if sepfigs, figure; end
    [axs rfweights{i} crweights{i} f{i} gof(i)] = allcones_plot(conerun, rasterrun, localmap, padfactor, 'ar', ar);        
    drawnow
    
end

ncones=[];
for i=1:length(rfweights)
    ncones(i)=length(rfweights{i});
end

save('/Volumes/Analysis/alex/peters_fits.mat', 'rfweights', 'crweights', 'f', 'gof', 'coneID', 'ncones', 'ncells')

figure
hold on
for i=1:length(rfweights)
    plot(rfweights{i}, f{i}.a.*crweights{i}, '.k');     
end
plot([-.2 1.2], [-0.2 1.2], 'k');
plot([-0.2 1.2], [0 0], 'k')
plot([0 0], [-0.2 1.2],'k')
axis([-0.2 1.2 -0.2 1.2])
xlabel('White Noise estimate')
ylabel('Contrast response estimate')



%% all cells binned cone input

% incl - only care for frames with reference cone at desired value
% excl - other cones should be above THR (OFF cells want negative  input)
% only - other cones should EACH have ABS value below THR
% str - other cones SUM should have ABS value below THR
% rand - assign random values, for control

load('/Volumes/Analysis/alex/peters_fits.mat', 'coneID', 'ncones')

mode='incl'; % 'incl', 'excl', 'rand', 'str', 'onl'
ninclude=0; % n of trials to include
plotit=false;
plotbybin=false;
plotalt=true;
plot_each_cone=false;
thr=0.1; % for "rest" in str, only, excl
bkgr=true;
alt_all=cell(1,length(ncells));
alt_std_all=alt_all;

% alt_coneID=cell(1, length(ncells));
for i = 1:length(ncells)
    
    conerun = eval(coneruns{ncells(i)});
    
    cellID=find(conerun.cell_ids==rgcs(i));
    
%     [val, ic]=sort(conerun.cones.weights(:,cellID), 'descend');
%     alt_coneID{i}=ic(1:ncones(i));
    
    myinputs=conerun.cone_inputs(:, coneID{i});
    myrate=double(conerun.spike_rate(cellID,:));
    if bkgr
        myrate=myrate-mean(myrate);
    end    

    mycone=zeros(10, 26, ncones(i));
    myconeSTD=mycone;
    alt=zeros(size(mycone,1), ncones(i));
    alt_std=alt;
    nentries=zeros(ncones(i),size(mycone,1));
    exclCones=cell(ncones(i), 10);
    
    if plotit
        figure
        set(gcf, 'Name', [dates(ncells(i), :), '  Cell ', int2str(rgcs(i)), '    ', mode])
        set(gcf, 'position', [1           1        1920        1105]);
    end
    
    for k=1:ncones(i)
        refcone=k;
        excones=setdiff(1:ncones(i),k);
        cnt=1;
        for bin=-0.5:0.1:0.4
            
            switch mode
                case 'incl'
                    tmp=find(myinputs(:,refcone)>= bin & myinputs(:,refcone)< (bin + 0.1));
                case 'excl'
                    tmp=find(myinputs(:,refcone)>= bin & myinputs(:,refcone)< (bin + 0.1) & ...
                        ~sum(myinputs(:,excones)<thr, 2));
                case 'str'
                    rest=sum(myinputs(:,excones),2);
                    tmp=find(myinputs(:,refcone)>= bin & myinputs(:,refcone)< (bin + 0.1) & ...
                        abs(rest)<abs(thr));
                case 'only'
                    tmp=find(myinputs(:,refcone)>= bin & myinputs(:,refcone)< (bin + 0.1) & ...
                        ~sum(abs(myinputs(:,excones))>thr, 2));
                case 'rand'
                    tmp=ceil(rand(100,1)*length(myrate-100));
            end
            
            tt=[];
            if ninclude
                trials=ninclude;
            else
                trials=length(tmp);
            end
            for j=1:trials
                if tmp(j)>6 && tmp(j)<length(myrate)-21
                    tt=[tt; myrate(tmp(j)-5:tmp(j)+20)];
                end
            end
            

            
            nentries(k,cnt)=size(tt,1);
            mycone(cnt,:,k)=mean(tt);
            alt_std(cnt,k)=std(tt(:,7));
%             mycone(cnt,:,k)=sum(tt);            
            alt(cnt,k)=mycone(cnt,7, k);
            
            
            cnt=cnt+1;
        end
        
        if plotit
            subplot(3,ceil(ncones(i)/3),k)
            plot(mycone(:,:,k)')
            title(int2str(nentries(k,:)))
            line([6 6], [-0.1 1.5], 'color', 'k')
            axis tight
            
        end
    end
    
    if plotbybin
        
        figure
        set(gcf, 'Name', [dates(ncells(i), :), '  Cell ', int2str(rgcs(i)) , ' BY BIN   ', mode])
        set(gcf, 'position', [1           1        1920        1105]);
        for j=1:10
            subplot(2,5, j)
            wr=squeeze(mycone(j,:,:));
            plot(wr)
            axis tight
            title(num2str(mean(nentries(:,j))))
        end
    end
    
    if plot_each_cone
        bin=-0.5:0.1:0.5;
        
        for k=1:ncones(i)
            figure
            set(gcf, 'Name', [dates(ncells(i), :), ', CONE ', int2str(k), ',  Cell ', int2str(rgcs(i)) , ' BY BIN   ', mode])
            set(gcf, 'position', [665   364   644   471]);
            
            mmax=max(mycone(:));
            mmin=min(mycone(:));
            for j=1:10
                subplot(2,5, j)
                wr=squeeze(mycone(j,:,k));
                plot(wr, 'r', 'linewidth', 3)
                axis([1 26 mmin mmax])
                line([1 26], [0 0], 'color', 'k')
                title(['[',num2str(bin(j)), ', ', num2str(bin(j)+0.1), '], n=', num2str(nentries(k,j))])
            end
        end
        
        tmp=zeros(ncones(i),26);
        for k=1:ncones(i)
            wr=squeeze(mycone(1:4,:,k));
            
            tmp(k,:)=sum(wr)/sum(nentries(k,1:4));
        end
        figure
        hold on
        plot(tmp', 'linewidth',2)
        title('summed, summed,  mean')
    end
    
    if plotalt
        
        figure
        set(gcf, 'Name', [dates(ncells(i), :), '  Cell ', int2str(rgcs(i)) , ' ALT   ', mode, '  ', int2str(ncones(i)), ' cones'])
        plot(alt, '-x')
        line([1 10], [0,0],'color','k')
        bin=-0.45:0.1:0.45;
        set(gca, 'xticklabel', {num2str(bin')})
        xlabel('middle of the bin')
        axis tight
    end
    
    alt_all{i}=alt;
    alt_std_all{i}=alt_std;
    
end


for i=1:4%length(ncells)
    figure
    subplot(1,2,1)
    plot(alt_all{i}, '-x')
    subplot(1,2,2)
    plot(alt_std_all{i}, '-x')
end


save('/Volumes/Analysis/alex/cr_white_noise_str_mode_thr_0.2.mat', 'alt_all', 'ncones')
save('/Volumes/Analysis/alex/cr_white_noise_incl_mode.mat', 'alt_all', 'mycone',  'ncones')



load('/Volumes/Analysis/alex/cr_white_noise_incl_mode.mat', 'alt_all', 'ncones')
for i = 1:length(ncells)
    figure
    set(gcf, 'Name', [dates(ncells(i), :), '  Cell ', int2str(rgcs(i)) , ' ALT   ', mode, '  ', int2str(ncones(i)), ' cones'])
    plot(alt_all{i}, '-x')
    line([1 10], [0,0],'color','k')
    bin=-0.45:0.1:0.45;
    set(gca, 'xticklabel', {num2str(bin')})
    xlabel('middle of the bin')
    axis tight
    title([dates(ncells(i), :), '  Cell ', int2str(rgcs(i)) , ' ALT   ', mode, '  ', int2str(ncones(i)), ' cones'])
end

% plot cone inputs histogram
figure
for i=1:9
    subplot(3,3,i)
    tmp=double(myinputs(:,i));
    [a, b]=hist(tmp, 25);
    fits=fit(b',a','Gauss1');
    hist(tmp, 25);
    hold on
    x=-0.5:0.01:0.5;
    y=fits.a1*exp(-((x-fits.b1)/fits.c1).^2);
    plot(x,y,'r', 'linewidth',2)
end


% get cone scaling from alt
myscaling=cell(1,length(ncells));
for i=1:length(ncells)
    tmp=alt_all{i}; % responses to 10 contrasts derived from WN
    [r,c]=find(tmp==max(tmp(:)));
    a=tmp(:, c);
    myscaling{i}(c)=1;
    for j=setdiff(1:size(tmp,2),c)
        [f gof] = fit(tmp(:, j), a, fittype({'x'}));
        myscaling{i}(j)=1/f.a;
    end
end

% plot cone weights from peter's crweights and mine rfweights

load('/Volumes/Analysis/alex/peters_fits.mat', 'rfweights', 'crweights', 'f', 'gof', 'coneID')

figure
for i=1:length(ncells)
    myrfweights=myscaling{i}./max(myscaling{i});    
    [f gof] = fit(crweights{i}, myrfweights', fittype({'x'}));
    
    subplot(4,4,i)
    plot(myrfweights, f.a.*crweights{i}, '.k'); hold on; plot([0 1], [0 1], 'k');
    plot([-.2 1.2], [-0.2 1.2], 'k');
    plot([-0.2 1.2], [0 0], 'k')
    plot([0 0], [-0.2 1.2],'k')
    axis([-0.2 1.2 -0.2 1.2])
    title(num2str(gof.rsquare))

end


figure
for i=1:length(ncells)
    myrfweights=myscaling{i}./max(myscaling{i});
    ic=myrfweights>0.33;
    [f gof] = fit(crweights{i}(ic), myrfweights(ic)', fittype({'x'}));
    plot(myrfweights(ic), f.a.*crweights{i}(ic), '.k'); hold on; plot([0 1], [0 1], 'k');
end
plot([-.2 1.2], [-0.2 1.2], 'k');
plot([-0.2 1.2], [0 0], 'k')
plot([0 0], [-0.2 1.2],'k')
axis([-0.2 1.2 -0.2 1.2])

load('/Volumes/Analysis/alex/peters_fits.mat', 'rfweights', 'crweights', 'f', 'gof', 'coneID')
hold on
for i=1:length(ncells)
    plot(rfweights{i}, f{i}.a.*crweights{i}, '.r');
end


figure
for i=1:length(ncells)
    myrfweights=rfweights{i}./max(rfweights{i});
    ic=myrfweights>0.33;
    [f gof] = fit(crweights{i}(ic), myrfweights(ic), fittype({'x'}));
    plot(myrfweights(ic), f.a.*crweights{i}(ic), '.k'); hold on; plot([0 1], [0 1], 'k');
end
plot([-.2 1.2], [-0.2 1.2], 'k');
plot([-0.2 1.2], [0 0], 'k')
plot([0 0], [-0.2 1.2],'k')
axis([-0.2 1.2 -0.2 1.2])




figure
for i=1:length(ncells)
    subplot(4,4,i)
    hold on
    a=crweights{i}./max(crweights{i});
    ic=a>0.1;
    a=a(ic);
    b=myscaling{i}'./max(myscaling{i});
    b=b(ic);
    c=rfweights{i};
    c=c(ic);
    [a, ic]=sort(a);
    b=b(ic);
    c=c(ic);
    
    err_ab=sqrt(sum((a-b).*(a-b)));
    err_bc=sqrt(sum((b-c).*(b-c)));
    err_ac=sqrt(sum((a-c).*(a-c)));

    plot(a, '*')
    plot(b, '+r')
    plot(c, '>g')
    if i==16
        legend('cr weights peter', 'cr weights WN', 'rf weight peter')
    end
    line([0 ncones(i)],[0 0], 'color', 'k')
    axis([0 ncones(i) -1 1.2])
    title(['CR vs WN ', num2str(err_ab, 2), '  WN vs RF ', num2str(err_bc, 2), '  CR vs RF', num2str(err_ac, 2)])
end


%% play with cancelling

load('/Volumes/Analysis/alex/peters_fits.mat')
load('/Volumes/Analysis/alex/myRFfit_incl_mode.mat') % load weight estimation (by my procedure)
load('/Volumes/Analysis/alex/subunits_weights.mat', 'subUnitID','subsweights')

% myrfweights=crweights;

bkgr=true; % subtract background firing
signBorder=0.25;
bin=[-0.5, -signBorder];% choose cone input range to look for
percentTolerance=0.15;
maxcones=5;
leg=false;
savepath='/Users/alexth/Desktop/single_cones/SUBUNITScancelling_signBorder_0.25_tolerance_0.15/myrfweights/';
if ~isdir(savepath)
    mkdir(savepath)
end


colors='cmgr';
maxcones=10
for i = 14:length(ncells) % loop through cells
    
    conerun = eval(coneruns{ncells(i)});    
    cellID=find(conerun.cell_ids==rgcs(i));
    run=conerun.names.rrs_prefix(end-6:end);
    
    
%     myData=dir(['/Volumes/Analysis/', dates(ncells(i),:), '/subunits/data*']);
    path2data=['/Volumes/Analysis/', dates(ncells(i),:), '/subunits/',run,'/anal/subunit/'];
    myData=dir([path2data, int2str(rgcs(i)), '*greedy*']);
    load([path2data, myData.name])
    
    locs=res.fit_SUB.locs_c;
    all_locs=conerun.cones.centers;
    myCone=coneID{i};
    refconessubunits=zeros(1, length(myCone));
    clear subunitCone
    for j=1:size(locs,1)
        subunitCone(j)=find(all_locs(:,1)==locs(j,1) & all_locs(:,2)==locs(j,2));
        tmp=find(myCone == subunitCone(j));
        if ~isempty(tmp)
            subUnitID{i}(j)=tmp;
        else
            subUnitID{i}(j)=0;
        end
    end
    
    presentInd=find(subUnitID{i});
    refconessubunits(subUnitID{i}(presentInd))=1;
    linInd=[];
    if size(res.fit_SUB.I_sc,1)~=size(res.fit_SUB.I_sc,2)
        subunits=find(sum(res.fit_SUB.I_sc')>1);
        for jj=1:length(subunits)
            subsCones=find(res.fit_SUB.I_sc(subunits(jj),:));
            [~,ia,~]=intersect(myCone, subunitCone(subsCones)); % cones with index ia are in subunit!
            refconessubunits(ia)=refconessubunits(ia)+jj;
            linInd=[linInd ia'];
        end
        linInd=[linInd setdiff(subUnitID{i}(presentInd),linInd)];
    else
        linInd=subUnitID{i}(presentInd);
    end
    

    % load raw cone inputs
    myinputs=conerun.cone_inputs(:, coneID{i}); 
    
    % multiply cone inputs by cone weight to equalize them, plot hists 
    myrfweights{i}=myrfweights{i}/max(myrfweights{i});
    figure
    for refcone=1:ncones(i)        
        myinputs(:,refcone)=myinputs(:,refcone)*myrfweights{i}(refcone);
        subplot(4, ceil(ncones(i)/4), refcone)
        hist(myinputs(:,refcone), -0.5:0.025:0.5)
        axis([-.5 0.5 0 Inf])
        title(['Cone ', int2str(refcone)])
    end
    
%     saveas(gcf, ['/Users/alexth/Desktop/single_cones/cancelling/corrected_cone_inputs_distr/Cell_',int2str(i),'_',dates(ncells(i),:), '_', int2str(rgcs(i)),'.pdf'], 'pdf')
    close(gcf)

    
    % load spiking response, subtract mean firing rate (not spont act!)
    myrate=double(conerun.spike_rate(cellID,:)); 
    if bkgr
        myrate=myrate-mean(myrate);
    end
        
    
            
%     tmpp=sort(myrfweights{i}, 'descend');
%     if length(tmpp)>maxcones
%         cones2take=find(myrfweights{i}>=tmpp(maxcones));
%         
%     else
%         cones2take=[1:ncones(i)]';
%     end

    cones2take=linInd(1:min(length(linInd),maxcones));

    numr=length(cones2take);
    
    figure
    set(gcf, 'Name', ['Cell ', int2str(rgcs(i)), ', ', int2str(ncones(i)),' cones'])
    set(gcf,'position',[1           1        1920        1105])
    cnt=1;
    col='k';
    col=repmat(col,1,length(cones2take)^2);
    frames=zeros(1, length(cones2take)^2);
    
    clear ylims
    % loop through cones
    for refcone=cones2take
%         excones=setdiff(cones2take,refcone); % indices of all other cones
        
        excones=cones2take;
        excones(excones==refcone)=[];
        
        % find instances of weighted cone inputs for the REF CONE within specified
        % range (bin parameter)
        tmp=find(myinputs(:,refcone)>= bin(1) & myinputs(:,refcone)< bin(2));
                
        % throw away instances too close to the beginning or to the end
        tmp(tmp<6)=[];tmp(tmp>length(myrate)-20)=[];
        
    

        for complCone=excones
                       
            refInputs=myinputs(:,refcone);
            otherInputs=myinputs(:,complCone);
            tmpSameAbs=find(abs(refInputs)<abs(otherInputs)*(1+percentTolerance) & abs(refInputs)>abs(otherInputs)*(1-percentTolerance));
            tmpbothPos=find(refInputs>signBorder & otherInputs>signBorder);
            tmpbothNeg=find(refInputs<-signBorder & otherInputs<-signBorder);
            tmpOpposite=find(refInputs<-signBorder & otherInputs>signBorder);
            tmpAllNeg=find(refInputs>=bin(1) & refInputs<=bin(2));
            tmpAllPos=find(refInputs<=-bin(1) & refInputs>=-bin(2));
            
            subsetSamePos=intersect(tmpSameAbs,tmpbothPos);
            subsetSameNeg=intersect(tmpSameAbs,tmpbothNeg);
            subsetOpposite=intersect(tmpSameAbs, tmpOpposite);
            
            tmpAllNeg(tmpAllNeg<6 | tmpAllNeg>length(myrate)-20)=[];
            tmpAllPos(tmpAllPos<6 | tmpAllPos>length(myrate)-20)=[];
            subsetSamePos(subsetSamePos<6 | subsetSamePos>length(myrate)-20)=[];
            subsetSameNeg(subsetSameNeg<6 | subsetSameNeg>length(myrate)-20)=[];
            subsetOpposite(subsetOpposite<6 | subsetOpposite>length(myrate)-20)=[];
         
            
            responseSubset=zeros(5, 25);
            for timePoints=1:25
                responseSubset(1, timePoints)=mean(myrate(subsetSamePos-6+timePoints)./abs(refInputs(subsetSamePos))');
                responseSubset(2, timePoints)=mean(myrate(subsetSameNeg-6+timePoints)./abs(refInputs(subsetSameNeg))');
                responseSubset(3, timePoints)=mean(myrate(subsetOpposite-6+timePoints)./abs(refInputs(subsetOpposite))');
                responseSubset(4, timePoints)=mean(myrate(tmpAllNeg-6+timePoints)./abs(refInputs(tmpAllNeg))'); % all values
                responseSubset(5, timePoints)=mean(myrate(tmpAllPos-6+timePoints)./abs(refInputs(tmpAllPos))'); % all values

            end
            
            [r,c]=ind2sub([numr,numr],cnt);
            if r==c
                subplot(numr,numr,cnt)
                plot(responseSubset(4:5,:)', 'linewidth', 2)
                axis tight
                ylims(1:2,cnt)=get(gca,'yLim');
                title(['C', int2str(refcone), ', w=', num2str(myrfweights{i}(refcone),2)...
                    ', peak=', num2str(ylims(2,cnt),2)])
                cnt=cnt+1;
                [r,c]=ind2sub([numr,numr],cnt);
            end
            
            if r~=c && refconessubunits(refcone)>1 && (refconessubunits(refcone)==refconessubunits(complCone))
                frames(cnt)=true;
                col(cnt)=colors(refconessubunits(refcone));
            else
                frames(cnt)=false;
                
            end
            
            if  r~=c
                subplot(numr,numr,cnt)
                plot(responseSubset(1:3,:)', 'linewidth', 2)
                axis tight
                ylims(1:2,cnt)=get(gca,'yLim');
                axis([1 25 -0.5 0.5])
                title(['C', int2str(refcone), ' vs C', int2str(complCone)])
                cnt=cnt+1;
                [r,c]=ind2sub([numr,numr],cnt);
                a=length(subsetSamePos);b=length(subsetSameNeg);d=length(subsetOpposite);
                if leg
                    legend(int2str([a b d]'));
                end
            end
            

            
            if r==c
                subplot(numr,numr,cnt)
                plot(responseSubset(4:5,:)', 'linewidth', 2)
                axis tight
                ylims(1:2,cnt)=get(gca,'yLim');
                title(['C', int2str(refcone), ', w=', num2str(myrfweights{i}(refcone),2)...
                    ', peak=', num2str(ylims(2,cnt),2)])
                cnt=cnt+1;
                [r,c]=ind2sub([numr,numr],cnt);
            end
        end
            
    end
    
    ymin=min(ylims(:));
    ymax=max(ylims(:));
    for subpl=1:cnt-1
        subplot(numr,numr,subpl)
        axis([1 25 ymin ymax]);
        line([7, 7], [ymin,ymax], 'color', 'k')
        line([1, 25], [0,0], 'color', 'k')
        set(gca,'xtick',0,'xticklabel','')
        
        
        if frames(subpl)
            rectangle('Position',[1,ymin,25,ymax-ymin], 'edgecolor',col(subpl), 'linewidth',2)
        end
    end
    drawnow
    
    saveas(gcf, [savepath, 'Cell_',int2str(i),'_',dates(ncells(i),:), '_', int2str(rgcs(i)),'.pdf'], 'pdf')
    close(gcf)
end



top_cells=[2 3];
linear_cells=[1 4 5 13 15 18 21];



%% center of mass
load('/Volumes/Analysis/alex/peters_fits.mat', 'rfweights', 'crweights', 'f', 'gof', 'coneID')

centers=cell(1, length(ncells));
cone_distances=centers;
centralConesN=1;
figure
col='k';
for i = 1:length(ncells)    
    
    conerun = eval(coneruns{ncells(i)});
    centers{i}=conerun.cones.centers(coneID{i}, :);
  
%     w=myscaling{i}'./max(myscaling{i});
    w=crweights{i}./max(crweights{i});
%     w=rfweights{i}./max(rfweights{i});


    [~, ik]=sort(w, 'descend');
    ik=ik(1:centralConesN);

    %com    

    x=centers{i}(:,1);
    y=centers{i}(:,2);

    
    x_com=sum(w(ik).*x(ik))/sum(w(ik));
    y_com=sum(w(ik).*y(ik))/sum(w(ik));
    
    cone_distances{i}=sqrt((x-x_com).^2+(y-y_com).^2);
    [~, ic]=sort(cone_distances{i});
    
    subplot(4,4,i)
    hold on
    plot(cone_distances{i}(ic),w(ic), col)
    axis([0 Inf 0 1])

end



%% Peter's approach to calculate raster weights

for i=1:length(ncells)
    conerun = eval(coneruns{ncells(i)});
    rasterrun = eval(rasterruns{ncells(i)});
    rasterrun.map=eval(['rasterrun.map' coneruns{ncells(i)}]);
    conerun.rgcs = rgcs(i);
    rasterrun.rgcs = rasterrun.map(get_cell_indices(conerun, conerun.rgcs));
    badregions = badregionss{i};
   
    
    padfactor = padfactors(i);
    badregions = badregionss{i};
    ar = ars(i);
    
    localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
    localmap = localmask.*rasterrun.stimulus.mapims{1};
 
    if ~isempty(badregions)
        for r = badregions, localmap(localmap == r) = 0; end
    end
    
    regions = setdiff(unique(localmap), 0);
    
    [crs crsx rasterhistx rasterhisty] = calc_allcones_cr(rasterrun, rasterrun.rgcs{1}, rasterrun.triggers(1:2:end));
    figure
    plot(crs(regions,:)', '-*')

end


figure
[p resnorm residual] = normcdfxscalesimple(crs(regions,:), crsx(regions,:), 'plot', true, 'title', false);
crweights = p(1:end-2)';

% crs is MxN response matrix, M is # cones, N is # contrasts
% crsx is MxN contrast matrix, with contrast values (negative for OFF cells)

% for cell 1 (2011-12-13-2, cell 1321)
crsx =[

   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0
   -1.4400   -1.0800   -0.7200   -0.3600         0];

crs = [

    3.3333    3.1667    2.6667    3.1667    2.5000
    3.0000    2.6667    3.3333    2.3333    2.5000
    4.3333    2.6667    3.6667    3.1667    2.5000
    2.6667    2.0000    2.5000    3.0000    2.5000
    2.8333    2.5000    3.6207    3.6667    2.5000
    2.3333    3.3333    2.6667    2.1667    2.5000
    2.5000    3.8333    2.3333    4.1667    2.5000
    3.3333    3.8333    3.1667    3.5000    2.5000
    1.5000    3.1667    2.6667    3.0000    2.5000
    1.8333    2.1667    2.6667    3.5000    2.5000
    3.0000    3.5000    2.5000    3.6667    2.5000
    1.8333    2.8333    2.8333    3.3333    2.5000
    7.5000    5.1667    4.0000    3.6667    2.5000
    5.1667    5.6667    5.3333    2.1667    2.5000
    8.8333    6.0000    3.6667    4.3333    2.5000
    9.3333    6.8333    4.6667    4.3333    2.5000
   10.3333    7.5000    4.6667    3.8333    2.5000
   12.0000    9.0000    6.3333    4.8333    2.5000
   11.1667    9.3333    5.5000    3.1667    2.5000
   15.1667   12.5000    7.6667    4.3333    2.5000];

regions = sparse([
    13
    14
    15
    16
    17
    18
    19
    20]);
%% binned approach for single cell from scratch (only cone inputs are needed), cone ID not provided, mode INCL

% params
conerun=d02s;
rgc=2389;

nconesCell=10;
bkgr=true;

% calculate responses for binned cone input
cellID=find(conerun.cell_ids==rgc);
[val, ic]=sort(conerun.cones.weights(:,cellID), 'descend');
myconeID=ic(1:nconesCell);
myinputs=conerun.cone_inputs(:, myconeID);
myrate=double(conerun.spike_rate(cellID,:));
if bkgr
    myrate=myrate-mean(myrate);
end
alt=zeros(10, nconesCell);
for refcone=1:nconesCell
    excones=setdiff(1:nconesCell,refcone);
    cnt=1;
    for bin=-0.5:0.1:0.4        

        tmp=find(myinputs(:,refcone)>= bin & myinputs(:,refcone)< (bin + 0.1));
        
        tmp(tmp<6 | tmp>(length(myrate)-20))=[];
        tt=zeros(1,25);
        for timePoint=1:25
            tt(timePoint)=mean(myrate(tmp-6+timePoint));
        end
        alt(cnt,refcone)=tt(7);
        cnt=cnt+1;
    end    
end
      
figure
plot(alt, '-x')
line([1 10], [0,0],'color','k')
bin=-0.45:0.1:0.45;
set(gca, 'xticklabel', {num2str(bin')})
xlabel('middle of the bin')
axis tight

% calculate relative cone weights (binned approach)
mycrs=alt(1:5, :)';
mycrsx=repmat([-0.4:0.1:0],nconesCell,1);
mycrs(:,5)=max(mean(mycrs(:,5)), 0);

[p resnorm residual] = normcdfxscalesimple(mycrs, mycrsx, 'plot', false, 'title', false);
myrfweights = p(1:end-2)';
myrfweights=myrfweights/max(myrfweights);

%plot cones with radius of weight index
locs=conerun.cones.centers(myconeID,:);
figure
viscircles(locs,myrfweights)
hold on
text(locs(:,1),locs(:,2),int2str([1:nconesCell]'));
axis ij
axis([0 600 0 600])


% cancelation analysis
%parameters
bkgr=true; % subtract background firing
signBorder=0.25;
bin=[-0.5, -signBorder];% choose cone input range to look for
percentTolerance=0.15;
maxcones=10;
leg=false;

%calculation
myinputs=conerun.cone_inputs(:, myconeID);

figure
for refcone=1:nconesCell
    myinputs(:,refcone)=myinputs(:,refcone)*myrfweights(refcone);
    subplot(4, ceil(nconesCell/4), refcone)
    hist(myinputs(:,refcone), -0.5:0.025:0.5)
    axis([-.5 0.5 0 Inf])
    title(['Cone ', int2str(refcone)])
end

myrate=double(conerun.spike_rate(cellID,:));
if bkgr
    myrate=myrate-mean(myrate);
end

tmpp=sort(myrfweights, 'descend');
if length(tmpp)>maxcones
    cones2take=find(myrfweights>=tmpp(maxcones));
    
else
    cones2take=[1:nconesCell]';
end

numr=length(cones2take);

figure
set(gcf, 'Name', ['Cell ', int2str(rgc), ', ', int2str(nconesCell),' cones'])
set(gcf,'position',[1           1        1920        1105])
cnt=1;

clear ylims
% loop through cones
for refcone=cones2take'
    excones=setdiff(cones2take,refcone); % indices of all other cones
    
    % find instances of weighted cone inputs for the REF CONE within specified
    % range (bin parameter)
    tmp=find(myinputs(:,refcone)>= bin(1) & myinputs(:,refcone)< bin(2));
    
    % throw away instances too close to the beginning or to the end
    tmp(tmp<6)=[];tmp(tmp>length(myrate)-20)=[];
    
    
    
    for complCone=excones'
%         
%         figure;plot(refInputs(tmpAllNeg),otherInputs(tmpAllNeg),'*')
%         tt=sum(myinputs(tmpAllNeg,[2 4:end])');
%         figure
%         hist(tt)
        
        refInputs=myinputs(:,refcone);
        otherInputs=myinputs(:,complCone);
        tmpSameAbs=find(abs(refInputs)<abs(otherInputs)*(1+percentTolerance) & abs(refInputs)>abs(otherInputs)*(1-percentTolerance));
        tmpbothPos=find(refInputs>signBorder & otherInputs>signBorder);
        tmpbothNeg=find(refInputs<-signBorder & otherInputs<-signBorder);
        tmpOpposite=find(refInputs<-signBorder & otherInputs>signBorder);
        tmpAllNeg=find(refInputs>=bin(1) & refInputs<=bin(2));
        tmpAllPos=find(refInputs<=-bin(1) & refInputs>=-bin(2));
        
        subsetSamePos=intersect(tmpSameAbs,tmpbothPos);
        subsetSameNeg=intersect(tmpSameAbs,tmpbothNeg);
        subsetOpposite=intersect(tmpSameAbs, tmpOpposite);
        
        tmpAllNeg(tmpAllNeg<6 | tmpAllNeg>length(myrate)-20)=[];
        tmpAllPos(tmpAllPos<6 | tmpAllPos>length(myrate)-20)=[];
        subsetSamePos(subsetSamePos<6 | subsetSamePos>length(myrate)-20)=[];
        subsetSameNeg(subsetSameNeg<6 | subsetSameNeg>length(myrate)-20)=[];
        subsetOpposite(subsetOpposite<6 | subsetOpposite>length(myrate)-20)=[];
        
        
%         tt=[];
%         for i=1:length(subsetOpposite)
%             tt=[tt; myrate((subsetOpposite(i)-3):(subsetOpposite(i)+6))];
%         end
%         figure
%         plot(tt')
%         figure
%         hist(tt(:,5), -1:0.5:5)
        
        
        responseSubset=zeros(5, 25);
        for timePoints=1:25
            responseSubset(1, timePoints)=mean(myrate(subsetSamePos-6+timePoints)./abs(refInputs(subsetSamePos))');
            responseSubset(2, timePoints)=mean(myrate(subsetSameNeg-6+timePoints)./abs(refInputs(subsetSameNeg))');
            responseSubset(3, timePoints)=mean(myrate(subsetOpposite-6+timePoints)./abs(refInputs(subsetOpposite))');
            responseSubset(4, timePoints)=mean(myrate(tmpAllNeg-6+timePoints)./abs(refInputs(tmpAllNeg))'); % all values
            responseSubset(5, timePoints)=mean(myrate(tmpAllPos-6+timePoints)./abs(refInputs(tmpAllPos))'); % all values
            
        end
        
        [r,c]=ind2sub([numr,numr],cnt);
        if r==c
            subplot(numr,numr,cnt)
            plot(responseSubset(4:5,:)', 'linewidth', 2)
            axis tight
            ylims(1:2,cnt)=get(gca,'yLim');
            title(['C', int2str(refcone), ', w=', num2str(myrfweights(refcone),2)...
                ', peak=', num2str(ylims(2,cnt),2)])
            cnt=cnt+1;
            [r,c]=ind2sub([numr,numr],cnt);
        end
        if  r~=c
            subplot(numr,numr,cnt)
            plot(responseSubset(1:3,:)', 'linewidth', 2)
            axis tight
            ylims(1:2,cnt)=get(gca,'yLim');
            axis([1 25 -0.5 0.5])
            title(['C', int2str(refcone), ' vs C', int2str(complCone)])
            cnt=cnt+1;
            [r,c]=ind2sub([numr,numr],cnt);
            a=length(subsetSamePos);b=length(subsetSameNeg);d=length(subsetOpposite);
            if leg
                legend(int2str([a b d]'));
            end
        end
        if r==c
            subplot(numr,numr,cnt)
            plot(responseSubset(4:5,:)', 'linewidth', 2)
            axis tight
            ylims(1:2,cnt)=get(gca,'yLim');
            title(['C', int2str(refcone), ', w=', num2str(myrfweights(refcone),2)...
                ', peak=', num2str(ylims(2,cnt),2)])
            cnt=cnt+1;
            [r,c]=ind2sub([numr,numr],cnt);
        end
    end
    
end

ymin=min(ylims(:));
ymax=max(ylims(:));
for subpl=1:cnt-1
    subplot(numr,numr,subpl)
    axis([1 25 ymin ymax]);
    line([7, 7], [ymin,ymax], 'color', 'k')
    line([1, 25], [0,0], 'color', 'k')
    set(gca,'xtick',0,'xticklabel','')
end


%% plot alt and save

load('/Volumes/Analysis/alex/cr_white_noise_incl_mode.mat', 'alt_all', 'mycone',  'ncones')
savepath='/Users/alexth/Desktop/single_cones/alt/';
for i=1:length(ncells)
    figure
    plot(alt_all{i}, '-x')
    line([1 10], [0 0], 'color', 'k')
    title(['Cell_',int2str(i),'_',dates(ncells(i),:), '_', int2str(rgcs(i))], 'Interpreter', 'None')
    saveas(gcf, [savepath, 'Cell_',int2str(i),'_',dates(ncells(i),:), '_', int2str(rgcs(i)),'.pdf'], 'pdf')
    close(gcf)
end


%% peter vs binned 
load('/Volumes/Analysis/alex/peters_fits.mat', 'rfweights', 'crweights', 'f', 'gof', 'coneID')
load('/Volumes/Analysis/alex/cr_white_noise_incl_mode.mat', 'alt_all', 'mycone',  'ncones')
load('/Volumes/Analysis/alex/cr_white_noise_str_mode_thr_0.2.mat', 'alt_all', 'ncones')


figure
sepplots=true;
twoplots=false;
thresh=0.25;
subset=1:21;
% subset=setdiff(subset, [linear_cells]);
% subset=linear_cells;
clear acc_peter acc_my
cnt=1;
for i=subset%1:length(ncells)
    mycrs=alt_all{i}(1:5, :)';
    mycrsx=repmat([-0.4:0.1:0],ncones(i),1);
    mycrs(:,5)=max(mean(mycrs(:,5)), 0);
    
    if sepplots
        subplot(5,5,i)
    end
    [p resnorm residual] = normcdfxscalesimple(mycrs, mycrsx, 'plot', false, 'title', false);
    myrfweights = p(1:end-2)';
    myrfweights=myrfweights/max(myrfweights);
    
    rr=crweights{i};
    
    ic=rr>max(rr)*thresh;
    [myf mygof] = fit(rr(ic), myrfweights(ic), fittype({'x'}));
    if twoplots
        subplot(1,2,1)
        title('Binned')
        line([0 1.3], [0 1.3],'color', 'k');
    end
    plot(myrfweights(ic), myf.a.*rr(ic), '.r');
    hold on; 
    axis([0 1.3 0 1.3])
    acc_my(cnt)=mygof.sse;
    
    if twoplots
        subplot(1,2,2)
        title('Standard')
        line([0 1.3], [0 1.3],'color', 'k');
    end
    [f gof] = fit(rr(ic), rfweights{i}(ic), fittype({'x'}));
    plot(rfweights{i}(ic), f.a.*rr(ic), '.k');
    axis([0 1.3 0 1.3])
    hold on; 
    acc_peter(cnt)=gof.sse;
    if ~twoplots
        title(['Binned ', num2str(mygof.sse,2), ', Standard ', num2str(gof.sse,2)])
        line([0 1.3], [0 1.3],'color', 'k');
    end
    if sepplots && i==1
        legend('Binned', 'Standard', 'location', 'best')    
    end
    if sepplots 
        line([0 1.3], [0 1.3],'color', 'k');
    end
    cnt=cnt+1;
end


figure
bar([acc_peter; acc_my]')
legend('Standard', 'binned', 'location', 'best')
title('SSE for standard and binned fit, linear cells, threshold 0.2')

figure
bar([1 5],[mean(acc_peter), mean(acc_my)]')
axis([-5 12 0 0.2])
set(gca,'xticklabel', {'standard', 'binned'})
title('Mean SSE for standard and binned fit, threshold 0')



top_cells=[2 3];
linear_cells=[1 4 5 13 15 18 21];% 13 questionable
bottom_cells=[6 7 8 9 10 11 12 14 16 17 19 20];
clear a
a(1,1)=mean(acc_peter(top_cells));
a(1,2)=mean(acc_my(top_cells));
a(2,1)=mean(acc_peter(linear_cells));
a(2,2)=mean(acc_my(linear_cells));
a(3,1)=mean(acc_peter(bottom_cells));
a(3,2)=mean(acc_my(bottom_cells));
figure
bar(a)
set(gca,'xticklabel', {'top', 'linear', 'bottom'})
legend('Standard', 'Binned')
title('SSE by groups of cells')


myrfweights=cell(1, length(ncells));
for i=1:length(ncells)
    mycrs=alt_all{i}(1:5, :)';
    mycrsx=repmat([-0.4:0.1:0],ncones(i),1);
    mycrs(:,5)=max(mean(mycrs(:,5)), 0);    
    [p resnorm residual] = normcdfxscalesimple(mycrs, mycrsx, 'plot', false);
    myrfweights{i} = p(1:end-2)';
    myrfweights{i} = myrfweights{i}/max(myrfweights{i});    
end
save('/Volumes/Analysis/alex/myRFfit_incl_mode.mat', 'myrfweights')

%% Jeremy's subunits model



dates = ['2011-12-13-2'; '2012-04-13-1'; '2012-09-06-0'; '2012-09-21-2'; '2012-08-21-2'];
coneruns = {'d08s', 'd06s', 'd04s', 'd09', 'd01s'};
ncells=[1 1 1 1 1 1 2 3 3 3 3 3 4 4 4 4 5 5 5 5 5];
rgcs = [1321 1351 2251 3586 4576 5162 6136 5103 5746 6031 6317 7022 5086 6346 6406 7696 391 2266 2828 5447 7682];
for i=1:size(coneruns,2)
    load(['/Volumes/Analysis/', dates(i,:), '/singlecones/dataruns.mat'])
end

subUnitID=cell(1, length(ncells));
subsweights=cell(1, length(ncells));
nsubunits=zeros(1,length(ncells));
for i=1:length(ncells)
    myData=dir(['/Volumes/Analysis/', dates(ncells(i),:), '/subunits/data*']);
    path2data=['/Volumes/Analysis/', dates(ncells(i),:), '/subunits/',myData.name,'/anal/subunit/'];
    myData=dir([path2data, int2str(rgcs(i)), '*greedy*']);
    load([path2data, myData.name])
    
    conerun = eval(coneruns{ncells(i)});
        
    locs=res.fit_SUB.locs_c;
    all_locs=conerun.cones.centers;
    myCone=coneID{i};
    clear subunitCone
    for j=1:size(locs,1)
        subunitCone(j)=find(all_locs(:,1)==locs(j,1) & all_locs(:,2)==locs(j,2));
        tmp=find(myCone == subunitCone(j));
        if ~isempty(tmp)
            subUnitID{i}(j)=tmp;
        else
            subUnitID{i}(j)=0;
        end
    end
    cnt=1;
    for j=1:size(res.fit_SUB.I_sc,1)
        subCones=res.fit_SUB.A_sc(j,find(res.fit_SUB.I_sc(j,:)));
        for k=1:length(subCones)
            subsweights{i}(cnt)=res.fit_SUB.B_s(j)*subCones(k);
            cnt=cnt+1;
            if k>1
                nsubunits(i)=nsubunits(i)+1;
            end
        end
    end
   
end

save('/Volumes/Analysis/alex/subunits_weights.mat', 'subUnitID','subsweights', 'nsubunits')


load('/Volumes/Analysis/alex/peters_fits.mat', 'rfweights', 'crweights', 'f', 'gof', 'coneID')
load('/Volumes/Analysis/alex/cr_white_noise_incl_mode.mat', 'ncones')
load('/Volumes/Analysis/alex/myRFfit_incl_mode.mat', 'myrfweights')
load('/Volumes/Analysis/alex/subunits_weights.mat', 'subUnitID','subsweights')


top_cells=[2 3];
linear_cells=[1 4 5 13 15 18 21];% 13 questionable
bottom_cells=[6 7 8 9 10 11 12 14 16 17 19 20];
figure
clear acc_peter acc_my acc_subs
cnt=1;
subset=1:21;
% subset=setdiff(subset, [linear_cells]);
% subset=linear_cells;
for i=subset
    subs=subsweights{i}'; % subunit model weights
    subindices=subUnitID{i}; % sununit cones indices
    tmp=find(subindices==0);
    if ~isempty(tmp)
        subs(tmp)=[];
        subindices(tmp)=[];
    end
    
    myRF=myrfweights{i}(subindices); % my fit
    cr=crweights{i}(subindices); % flash response fit  
    peterRF=rfweights{i}(subindices); % peter's fit
    
    myRF=myRF/max(myRF); % my fit
    peterRF=peterRF/max(peterRF); % peter's fit
    subs = subs/max(subs);
    
    
    [myf mygof] = fit(cr, myRF, fittype({'x'})); % my fit vs flash
    [peterf petergof] = fit(cr, peterRF, fittype({'x'})); % my fit vs flash
    [subsf subsgof] = fit(cr, subs, fittype({'x'})); % my fit vs flash
    
    subplot(5,5,i)
    hold on;
    plot(myRF, myf.a.*cr, 'or');
    plot(peterRF, peterf.a.*cr, 'xk');
    plot(subs, subsf.a.*cr, '+b');
     
   
    acc_my(cnt)=mygof.sse;
    acc_peter(cnt)=petergof.sse;
    acc_subs(cnt)=subsgof.sse;
    
    title(['Binned ', num2str(mygof.sse,2), ', Standard ', num2str(petergof.sse,2), ', Subunits ', num2str(subsgof.sse,2)])

    if i==1
        legend('Binned', 'Standard','subunits', 'location', 'best')    
    end
    axis([0 1.3 0 1.3])
    line([0 1.3], [0 1.3],'color', 'k');
    
    cnt=cnt+1;
end


figure
bar([acc_peter; acc_my; acc_subs]')
legend('Standard', 'binned','subunits', 'location', 'best')
title('SSE for standard, binned, subunits fit, LINEAR CELLS')

figure
subplot(1,2,1)
bar([1 5 9],[mean(acc_peter), mean(acc_my), mean(acc_subs)]')
axis([-5 16 0 0.16])
set(gca,'xticklabel', {'peter', 'alex', 'jeremy'}, 'fontsize',16, 'fontweight', 'b')
title('Mean SSE, LINEAR cells', 'fontsize',16, 'fontweight', 'b')

subplot(1,2,2)
bar([1 5 9],[mean(acc_peter), mean(acc_my), mean(acc_subs)]')
axis([-5 16 0 0.18])
set(gca,'xticklabel', {'peter', 'alex', 'jeremy'}, 'fontsize',16, 'fontweight', 'b')
title('Mean SSE, NONLINEAR cells', 'fontsize',16, 'fontweight', 'b')



%% do Peter's style but different stuff

sepfigs = false;
if ~sepfigs, figure; end
clear rfweights crweights f gof

mode_weights = 'sum';
for i = 1:length(ncells)
    
    conerun = eval(coneruns{ncells(i)});
    rasterrun = eval(rasterruns{ncells(i)});
    rasterrun.map=eval(['rasterrun.map' coneruns{ncells(i)}]);
    conerun.rgcs = rgcs(i);
    
    rasterrun.rgcs = rasterrun.map(get_cell_indices(conerun, conerun.rgcs));
    
    padfactor = padfactors(i);
    badregions = badregionss{i};
    ar = ars(i);
   
    localmask = stimmasklocal(conerun, conerun.rgcs, rasterrun.stimulus.mapims{1}, 'az_pad_factor', padfactor);
    localmap = localmask.*rasterrun.stimulus.mapims{1};
    if ~isempty(badregions)
        for r = badregions, localmap(localmap == r) = 0; end
    end
    
    if sepfigs
        figure;
        set(gcf, 'position', [94 677 560 420])
    end
    [axs rfweights{i} crweights{i} f{i} gof(i)] = allcones_plot_play(conerun, rasterrun, localmap, padfactor, 'ar', ar, 'crweights', mode_weights);        
    drawnow
    
end


save('/Volumes/Analysis/alex/sum_resp_fits.mat', 'rfweights', 'crweights', 'f', 'gof')
save('/Volumes/Analysis/alex/max_resp_fits.mat', 'rfweights', 'crweights', 'f', 'gof')

% sum vs Peter (cr weights only)
load('/Volumes/Analysis/alex/peters_fits.mat', 'rfweights', 'crweights', 'f', 'gof')
cr_peter=crweights;
f_peter=f;
load('/Volumes/Analysis/alex/sum_resp_fits.mat', 'rfweights', 'crweights', 'f', 'gof')
figure
hold on
for i=1:length(crweights)
    
    subplot(5,5,i)
    plot(cr_peter{i}/max(cr_peter{i}), crweights{i}, '.k');
    hold on
    plot([-.2 1.3], [-0.2 1.3], 'k');
    plot([-0.2 1.3], [0 0], 'k')
    plot([0 0], [-0.2 1.3],'k')
    axis([-0.2 1.3 -0.2 1.3])
    xlabel('Standard')
    ylabel('Sum response')
    title('Standard vs sum')

end

figure
subplot(2,2,1)
hold on
for i=1:length(top_cells)
    plot(cr_peter{top_cells(i)}/max(cr_peter{top_cells(i)}), crweights{top_cells(i)}, '.k');
    plot([-.2 1.3], [-0.2 1.3], 'k');
    plot([-0.2 1.3], [0 0], 'k')
    plot([0 0], [-0.2 1.3],'k')
    axis([-0.2 1.3 -0.2 1.3])
    xlabel('Standard')
    ylabel('Sum response')
    title('Standard vs sum, top')
end
subplot(2,2,2)
hold on
for i=1:length(linear_cells)
    plot(cr_peter{linear_cells(i)}/max(cr_peter{linear_cells(i)}), crweights{linear_cells(i)}, '.k');
    plot([-.2 1.3], [-0.2 1.3], 'k');
    plot([-0.2 1.3], [0 0], 'k')
    plot([0 0], [-0.2 1.3],'k')
    axis([-0.2 1.3 -0.2 1.3])
    xlabel('Standard')
    ylabel('Sum response')
    title('Standard vs sum, linear')
end
subplot(2,2,3)
hold on
for i=1:length(bottom_cells)
    plot(cr_peter{bottom_cells(i)}/max(cr_peter{bottom_cells(i)}), crweights{bottom_cells(i)}, '.k');
    plot([-.2 1.3], [-0.2 1.3], 'k');
    plot([-0.2 1.3], [0 0], 'k')
    plot([0 0], [-0.2 1.3],'k')
    axis([-0.2 1.3 -0.2 1.3])
    xlabel('Standard')
    ylabel('Sum response')
    title('Standard vs sum, bottom')
end




load('/Volumes/Analysis/alex/sum_resp_fits.mat', 'rfweights', 'crweights', 'f', 'gof')
figure
hold on
for i=1:length(rfweights)
    plot(rfweights{i}, f{i}.a.*crweights{i}, '.k');  
    summed_sse(i)=gof(i).sse;
end
plot([-.2 1.2], [-0.2 1.2], 'k');
plot([-0.2 1.2], [0 0], 'k')
plot([0 0], [-0.2 1.2],'k')
axis([-0.2 1.2 -0.2 1.2])
title('sum response fit')





load('/Volumes/Analysis/alex/peters_fits.mat', 'rfweights', 'crweights', 'f', 'gof')


figure
hold on
for i=1:length(rfweights)
    tmp=rfweights{i}/max(rfweights{i});
    crweights{i}=crweights{i}(tmp>0);
    tmp=tmp(tmp>0);
    [f gof] = fit(crweights{i},  tmp, fittype({'x'}));
    subplot(4,4,i)
    hold on
    plot(tmp, f.a.*crweights{i}, '.k'); hold on; plot([0 1], [0 1], 'k');
    title(num2str(gof.sse,2))
    axis([0 1.5 0 1.5])
end






%%

figure
hold on
for i=1:length(crweights)
    tmp=cr_peter{i}/max(cr_peter{i});
    [f gof] = fit(tmp,  crweights{i}, fittype({'x'}));
%     subplot(4,4,i)
    plot(crweights{i}, f.a.*tmp, '.k'); hold on; plot([0 1], [0 1], 'k');

%     plot(cr_peter{i}/max(cr_peter{i}), crweights{i}, '.k');
end
plot([-.2 1.3], [-0.2 1.3], 'k');
plot([-0.2 1.3], [0 0], 'k')
plot([0 0], [-0.2 1.3],'k')
axis([-0.2 1.3 -0.2 1.3])
xlabel('Peter')
ylabel('Sum response')
title('peter vs sum')


load('/Volumes/Analysis/alex/peters_fits.mat', 'rfweights', 'crweights', 'f', 'gof')
aa=[];bb=[];
for i=1:length(crweights)
    aa=[aa; sort(crweights{i}/max(crweights{i}))];
    bb=[bb; sort(rfweights{i}/max(rfweights{i}))];
end
figure
hist(aa, 0.05:0.1:1)
figure
hist(bb, 0.05:0.1:1)
figure
plot(sort(aa))
hold on
plot(sort(bb), 'r')

x = 0:0.01:2;
y = cdf('Norm',x,1,0.5);
figure
plot(-x, y)

hold on
y=0.5*(0:0.01:2);
plot(-x, y, 'r')



y = cumsum(cdf('Norm',0:0.01:2,1,0.5));
figure
plot(y)



