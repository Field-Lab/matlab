datarun = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data008-from-d08_11/data008-from-d08_11');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = load_neurons(datarun);
datarun = set_polarities(datarun);

wn_movie_name = 'BW-2-6-0.48-11111-300x300-60.35.xml';
[inputs, refresh, duration] = get_wn_movie_ath(datarun, wn_movie_name);

vormap = load('/Volumes/Data/2011-12-13-2/Visual/2011-12-13-2_f04_vorcones/map-0000.txt');

vorrun = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data009-from-d08_11/data009-from-d08_11');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);

[inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-937x1-60.35.xml');


sta_params.length = 3;
sta_params.offset = 0;
offset = 0;
sta_length=2;



for datarunID = 1:311
    
    visionID = datarun.cell_ids(datarunID);
    [folder, my_type] = find_cell_type(datarun, visionID);
    
    if ~strcmp(folder,'crap') && ~strcmp(folder,'duplicates')
        
        % full voronoi
        my_sta=zeros(size(inputs_v,1),sta_params.length);
        spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
        spike_rate=zeros(size(inputs_v,2),1);
        spikes_tmp = spikes(spikes>sta_params.length-sta_params.offset);
        full_voronoi_sta = zeros(600,600);
        nspikes = length(spikes_tmp);
        while ~isempty(spikes_tmp)
            [~, ia, ~] = unique(spikes_tmp);
            spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
            for j=1:sta_params.length
                my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                    sum(inputs_v(:,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
            end
            spikes_tmp(ia)=[];
        end
        % my_sta=my_sta/nspikes;
        
        % get voronoi map sta
        tmp_map = vormap;
        tt=0:max(tmp_map(:));
        for i=1:937
            [a, b] = find(tmp_map==tt(i+1));
            if ~isempty(a)
                for j = 1:length(a)
                    full_voronoi_sta(a(j),b(j)) = my_sta(i,3);
                end
            end
        end
        
        
        val = max(abs(my_sta(:,3)));
        ind = find(abs(full_voronoi_sta)==val, 1);
        [a,b] = ind2sub([600,600], ind);
        a = round(a/2);
        b = round(b/2);
        
        border = 50;
        myinds = [];
        my_as = a-border:a+border;
        my_as(my_as<1 | my_as>300) = [];
        my_bs = b-30:b+border;
        my_bs(my_bs<1 | my_bs>300) = [];
        for i = my_as
            for j=my_bs
                myinds = [myinds sub2ind([300, 300],i, j)];
            end
        end
        
        % do single cone WN
        
        spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/refresh); % spikes in frames
        
        spike_rate=zeros(duration,1);
        my_sta=zeros(length(myinds),sta_length);
        ksta = [sta_params.length-sta_params.offset, 100:100:18000];
        sta_tmp = zeros(border*2+1,border*2+1,length(ksta)-1);
        
        for cnt1 = 1:length(ksta)-1            
            spike_tmp = spikes(spikes>ksta(cnt1) &spikes<=ksta(cnt1+1));
            spike_tmp(spike_tmp<sta_length-offset)=[];
            nsp(cnt1) = length(spikes_tmp);
            while ~isempty(spike_tmp)
                [c, ia, ic] = unique(spike_tmp);
                spike_rate(spike_tmp(ia))=spike_rate(spike_tmp(ia))+1;
                for j=1:sta_length
                    my_sta(:,sta_length-j+1)=my_sta(:,sta_length-j+1)+sum(inputs(myinds,spike_tmp(ia)-sta_length+j+offset),2);
                end
                spike_tmp(ia)=[];
            end
            tmp = reshape(my_sta,length(my_bs),length(my_as),2);
            
            sta_tmp(:,:,cnt1) = tmp(:,:,2)';
        end
        
        save(['/Volumes/Analysis/2011-12-13-2/sta_time/single_cone_',int2str(datarunID)],'sta_tmp', 'nsp');
    end
end
