datarun = load_data('/Volumes/Analysis/2016-01-05-1/d06-21-norefit/data014/data014');
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

plot(diff(datarun.triggers(1:10)))



starun = load_data('/Volumes/Analysis/2016-01-05-1/d06-21-norefit/data013/data013');
starun = load_params(starun,'verbose',1);
starun = load_neurons(starun);



a = cell(1, length(variable_parameters));
stimpol = zeros(1, length(variable_parameters));
for i=1:length(variable_parameters)
    if variable_parameters(i).rgb(1) == -0.48
        stimpol(i) = -1;
    else
        stimpol(i) = 1;
    end
    a{i} = variable_parameters(i).index_map;
end

unique_stims = unique(a);
tot_stim = zeros(length(unique_stims),10);
for i=1:length(variable_parameters)
    t = find(strcmp(unique_stims,a{i}));
    p = find(tot_stim(t,:)~=0);
    if isempty(p)
        tot_stim(t,1) = i*stimpol(i);
    else
        tot_stim(t,p(end)+1) = i*stimpol(i);
    end
end

cellID = 53;
myspikes = cell(1,length(variable_parameters));
cnt = 1;
for i=1:2:length(variable_parameters)*2-2
    myspikes{cnt} = datarun.spikes{cellID}(datarun.spikes{cellID}>=datarun.triggers(i)-1 & datarun.spikes{cellID}<datarun.triggers(i+2));
    myspikes{cnt} = myspikes{cnt} - datarun.triggers(i);
    cnt=cnt+1;
end

figure
hold on
cnt=1;
for i=1:1500
    if strcmp(a{i},'on_parasol_61_sd_3');
        if stimpol(i)>0
            plot(myspikes{i}, ones(length(myspikes{i}))+cnt, 'xb')
        else
            plot(myspikes{i}, ones(length(myspikes{i}))+cnt, 'xk')
        end
        cnt=cnt+1;
    end
end

figure
plot(datarun.spikes{cellID}(1:2000))

figure
cc=1;
for j=91:100
    tmp = tot_stim(j,:);
    a = find(tmp>0);
    subplot(3,4,cc)
    hold on
    cnt = 1;
    for i=a
        plot(myspikes{abs(tot_stim(j,i))}, ones(1, length(myspikes{abs(tot_stim(j,i))}))*cnt, '.b')
        cnt = cnt+1;
    end
    axis([-1.5 5.5 0 6])
    title(unique_stims{j}, 'interpreter', 'none')
    cc=cc+1;
end

figure
plot(diff(datarun.triggers))


