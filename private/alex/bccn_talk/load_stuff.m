function [inputs, inputs_wn, inputs_nsem, meta, spikes, sta] = load_stuff(date, ndf, my_movie)

if strcmp(date, '2015-03-09-2')
    path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/';   
elseif strcmp(date, '2015-08-17-1')
    path2data = '/Volumes/Analysis/2015-08-17-1/d01-29-norefit/';
end

[nsem_data, sta_data, wn_data, movie_name]=choose_data(date, ndf);
nsemrun = load_other_data(path2data, nsem_data);
wnrun = load_other_data(path2data, wn_data);
starun = load_sta_data(path2data, sta_data);
[inputs, refresh, ~] = get_wn_movie_ath(starun, movie_name);
inputs = logical((inputs+0.48)/0.96);
[inputs_wn, ~, ~] = get_wn_movie_ath(wnrun, movie_name);
inputs_wn = logical((inputs_wn+0.48)/0.96);

stix_size=320/sqrt(size(inputs,1));
inputs_nsem = scale_movie(my_movie, stix_size);

meta.triggers.wn = wnrun.triggers;
meta.triggers.nsem = nsemrun.triggers;
meta.refresh = refresh;

sta = zeros(size(inputs,1), length(starun.cell_ids));

for i=1:length(starun.cell_ids)
    spikes.nsem{i} = nsemrun.spikes{i};
    spikes.wn{i} = wnrun.spikes{i};
    tmp = starun.spikes{i} - starun.triggers(1);
    tmp(tmp<0) = [];
    spikes.sta{i} = tmp;    
    pol = starun.stas.polarities{i};
    if ~isempty(pol) && pol~=0        
        sta_tmp = pol*squeeze(starun.stas.stas{i});
        sta_tmp = reshape(sta_tmp, [],30);
        sta_tmp = sta_tmp(:,11:end);
        [~,max_loc] = find(sta_tmp==max(sta_tmp(:)));
        a = robust_std(sta_tmp(:,max_loc(1)));
        mean_tc = mean(sta_tmp(sta_tmp(:,max_loc(1))>a*4,:),1);
        sta(:,i) = sum(sta_tmp.*repmat(mean_tc,size(sta_tmp,1),1),2);
    end
end


    function [nsem_data, sta_data, wn_data, movie_name]=choose_data(date,ndf)
        if strcmp(date, '2015-03-09-2')
            switch ndf
                case 0
                    nsem_data = 'data024';
                    sta_data = 'data026';
                    wn_data = 'data025';
                    movie_name = 'BW-4-2-0.48-11111-80x80.xml';
                case 1
                    nsem_data = 'data020';
                    sta_data = 'data022';
                    wn_data = 'data021';
                    movie_name = 'BW-8-6-0.48-11111-40x40.xml';
                case 2
                    nsem_data = 'data016';
                    sta_data = 'data018';
                    wn_data = 'data017';
                    movie_name = 'BW-8-6-0.48-11111-40x40.xml';
                    
                case 3
                    nsem_data = 'data012';
                    sta_data = 'data014';
                    wn_data = 'data013';
                    movie_name = 'BW-8-8-0.48-11111-40x40.xml';
                case 4
                    nsem_data = 'data009';
                    sta_data = 'data008';
                    wn_data = 'data010';
                    movie_name = 'BW-10-8-0.48-11111-32x32.xml';
            end
            
        elseif strcmp(date, '2015-08-17-1')
            switch ndf
                case 2
                    nsem_data = 'data014';
                    sta_data = 'data011';
                    wn_data = 'data012';
                    movie_name = 'BW-10-6-0.48-11111-32x32.xml';
                case 3
                    nsem_data = 'data010';
                    sta_data = 'data007';
                    wn_data = 'data008';
                    movie_name = 'BW-16-8-0.48-11111-20x20.xml';
                case 4
                    nsem_data = 'data005';
                    sta_data = 'data006';
                    wn_data = 'data004';
                    movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            end
        end
    end

    function otherrun = load_other_data(path2data, my_data)
        otherrun = load_data(fullfile(path2data, my_data, my_data));
        otherrun = load_params(otherrun,'verbose',1);
        otherrun = load_neurons(otherrun);
    end

    function starun = load_sta_data(path2data, my_data)
        starun = load_data(fullfile(path2data, my_data, my_data));
        starun = load_params(starun,'verbose',1);
        starun = load_neurons(starun);
        starun = load_sta(starun, 'load_sta', 'all');
        starun = set_polarities(starun);
    end

end
