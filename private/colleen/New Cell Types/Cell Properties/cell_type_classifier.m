% clustering cell types
clear
close all
%2007-05-01-2 is dark
%2007-08-21-4 has bad OFF large
%2009-04-13-1 has super large OFF 1
%2008-04-08-2 had ON large 2 that don't cross 0
%    '2005-04-14-0'     'data002';; outlier
   % '2008-04-22-5';    'data003'; outlier
   
     %  '2009-04-13-7';     'data000'; outlier
%          '2008-08-27-6';    'data009';

   % '2007-08-21-4';     'data000'; only 1 type

%     '2008-04-08-2';     'data003'; only 1 type

%     '2009-12-03-2';     'data000';
%     '2010-11-22-4';     'data000-nwpca';

%     '2013-05-28-9';      'data000';

[~, txt] = xlsread('/Users/colleen/Documents/Large Cell Data ARVO.xlsx');
    data= [];

for t= 1:size(txt,1)
    % piece = txt(j,1:3);
    run_opts.date{t}=strtrim(txt{t,1}); % one slash at the end
    temp1 = strtrim(txt{t,2});
    temp2 =  strtrim(txt{t,3});
    run_opts.concatname=temp1; % Name (or modified name) of run, no slashes\
    run_opts.dataname{t} = temp2;
    
    % Sometimes the data has two versions of the concate name
    run_opts.file_name{t} = [run_opts.date{t}, '/', run_opts.concatname, '/',  run_opts.dataname{t}];
    % run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];
    
    
    run_opts.save_location_root = '/Volumes/Lab/Users/crhoades/Cell Properties/';
    % Where to save
    run_opts.filepath{t}= [run_opts.save_location_root, run_opts.date{t}, '/', run_opts.concatname, '/', run_opts.dataname{t}, '/'];
    


% 
% date = {
%     };
% 
% concate = {
%     };
episolon = 0.0001;
% for t = 1:size(date,1)-4
    clear match
    filename = ['/Volumes/Lab/Users/crhoades/Cell Properties/', run_opts.file_name{t}];
    classes= {'OFF Parasol', 'ON Parasol', 'ON Midget', 'OFF Midget', 'ON Large 1', 'ON Large 2', 'ON Large 3', 'ON Large 4', 'OFF Large 1', 'OFF Large 2', 'OFF Large 3', 'OFF Large 4'};
    classes= {'OFF Parasol', 'ON Parasol',  'ON Large 1', 'ON Large 2', 'ON Large 3', 'ON Large 4', 'OFF Large 1', 'OFF Large 2', 'OFF Large 3', 'OFF Large 4'};
%     classes = {'ON Parasol', 'ON large 1', 'ON large 2'};
    classes = {'ON Parasol', 'ON large 1', 'ON large 2', 'OFF parasol', 'OFF large 1', 'OFF large 2'};

    D = dir(filename);
    colors = [];
    col = hsv(size(classes,2));
    counter = 1;
    for c = 1:size(classes,2)
        for i =1:size(D,1)
            match(c,i) = strcmpi(D(i).name, classes{c});
        end
        if sum(match(c,:))>0
            
            file_path{c} = [filename, '/', D(match(c,:)).name, '/', 'output.mat'];
            output = load(file_path{c});
            try
            new_data = [c*ones(size(cell2mat(struct2cell(output.output.parameters)'),1),1), cell2mat(struct2cell(output.output.parameters)')];
            catch
                output.output.parameters.pca1 = nan(length(output.output.parameters.fr));
                output.output.parameters.pca2 = nan(length(output.output.parameters.fr));
                new_data = [c*ones(size(cell2mat(struct2cell(output.output.parameters)'),1),1), cell2mat(struct2cell(output.output.parameters)')];
            end
            new_data = [new_data]; % remove acf
            if c == 1
                ref_on_parasol = new_data;
                    new_data(:,7) = 1;
                new_data = [repmat(t, size(new_data,1),1), new_data];

                data = [data; new_data(:, 1:9)];
            elseif c == 4
                ref_off_parasol = new_data;
                    new_data(:,7) = -1;
                new_data = [repmat(t, size(new_data,1),1), new_data];

                data = [data; new_data(:, 1:9)];
            end
            
            if c == 2 || c == 3 
                real_ind = boolean(floor(sum(~isnan(ref_on_parasol),2)/size(ref_on_parasol,2)));
                new_data(:,2:8) = new_data(:,2:8)./repmat(median(ref_on_parasol(real_ind, 2:8)), size(new_data,1),1);
                mean_amp = median(ref_on_parasol(real_ind, 8));
                if mean_amp < 0 
                    new_data(:,8) = -1;
                else
                    new_data(:,8) = 1;
                end
                new_data = [repmat(t, size(new_data,1),1), new_data];

                data = [data; new_data(:, 1:9)];

            elseif c == 5|| c == 6 
                real_ind = boolean(floor(sum(~isnan(ref_off_parasol),2)/size(ref_off_parasol,2)));
                new_data(:,2:8) = new_data(:,2:8)./repmat(median(ref_off_parasol(real_ind, 2:8)), size(new_data,1),1);
                 mean_amp = median(ref_off_parasol(real_ind, 8));
                if mean_amp < 0 
                   new_data(:,8) = -1;
                else
                    new_data(:,8) = 1;
                end
                new_data = [repmat(t, size(new_data,1),1), new_data];
                data = [data; new_data(:, 1:9)];

            end
            
            %         new_data_unit = new_data/std(new_data);
%             data = [data; new_data];
    %         colors = [colors ;repmat(col(counter, :), size(cell2mat(struct2cell(output.output.parameters)'),1),1)];

            %         colors = [colors ;repmat(col(counter, :), size(cell2mat(struct2cell(output.output.parameters)'),1),1)];
            counter = counter+1;
        end

    end
   
end


data_large = data(data(:,2) == 2 | data(:,2) == 3 | data(:,2) == 5 | data(:,2) == 6 | data(:,2) == 7| data(:,2) == 8 | data(:,2) == 9, :);
data_rows = boolean(floor(sum(~isnan(data_large),2)/size(data_large,2)));
% data_large(:, 3:end) = data_large(:,3:end)./repmat(max(abs(data_large(data_rows,3:end))), size(data_large,1),1);
% norm_data = data_large(:,3:end)./(repmat(std(data_large(data_rows,3:end)), size(data_large,1),1)+episolon);
% data_large(data_large(:,18)<-5, 18) = -5;
% data_large(data_large(:,17)>5, 17) = 5;

norm_data = data_large;%[data_large(:,1:2), norm_data];


[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_data(:,3:end), 'NumComponents',3);


% xdata = data_large(:,3:4);
% group = data_large(:,2);
% svmStruct = svmtrain(xdata, group, 'ShowPlot', true)
% [b,dev,stats] = glmfit(data_large(:,3:4),data_large(:,2)-2,'binomial','link','identity')
% % [idx, z] = mnrfit(norm_data(:,3:end), norm_data(:,2));
% y_fit = glmval(b,data_large(:,3:4), 'identity');
% figure
% plot(data_large(:,3:4), data_large(:,2)-2, 'o', data_large(:,3:4)-2, y_fit, 'o')
test = SCORE;
% idx = kmeans(norm_data(:,2:end),length(unique(norm_data(:,1))));
% idx = kmeans(norm_data,length(unique(norm_data(:,1))));
% idx = kmeans(test,4);

figure
scatter3(test(:,1), test(:,2),test(:,3));
title('PC 1,2, and 3')

figure;
% subplot(3,1,1)
set(gcf, 'Position', [10 10 1120 840])

h = gscatter(test(:,1), test(:,2),norm_data(:,1:2));

S = unique(norm_data(:,1:2), 'rows');
unique_dates = unique(norm_data(:,1));

unique_classes = unique(norm_data(:,2));
markers = {'+','o','*','square', 'hexagram','diamond','^','v','>','<','pentagram','x','+','o','*','square', 'hexagram','diamond','^','v','>','<','pentagram','x'};
% colors = hsv(size(unique(S(:,1)),1));
colors = hsv(length(unique_classes));

l= findobj(gcf,'tag','legend'); 
for n = 1:length(h)
   set(h(n), 'Marker', markers{unique_dates == S(n,1)})
      set(h(n), 'MarkerSize', 9);

   set(h(n), 'MarkerFaceColor', colors(unique_classes == S(n,2), :))
      set(h(n), 'MarkerEdgeColor', colors(unique_classes == S(n,2), :))
LegendString{n} = [[run_opts.date{S(n,1)}, ' ', run_opts.dataname{S(n,1)}(1:7)] ' , ' classes{S(n,2)}];
end
set(l, 'String', LegendString)
title('PC 1 and 2')
xlabel('PC 1')
ylabel('PC 2')
xaxis = get(gca, 'xlim');
set(gca, 'xlim', [xaxis(1), xaxis(2)+2])

set(l,'location','southeast');


figure;
set(gcf, 'Position', [10 10 1120 840])
% subplot(3,1,1)
h= gscatter(test(:,1), test(:,3),norm_data(:,1:2));

l= findobj(gcf,'tag','legend'); 
for n = 1:length(h)
   set(h(n), 'Marker', markers{unique_dates == S(n,1)})
         set(h(n), 'MarkerSize', 9);
   set(h(n), 'MarkerFaceColor', colors(unique_classes == S(n,2), :))
      set(h(n), 'MarkerEdgeColor', colors(unique_classes == S(n,2), :))
LegendString{n} = [run_opts.date{S(n,1)} ' , ' classes{S(n,2)}];
end
set(l, 'String', LegendString)
title('PC 1 and 3')
xlabel('PC 1')
ylabel('PC 3')
xaxis = get(gca, 'xlim');
set(gca, 'xlim', [xaxis(1), xaxis(2)+2])
set(l,'location','southeast');

set(l,'location','southeast');


figure;
set(gcf, 'Position', [10 10 1120 840])

% subplot(3,1,1)
h = gscatter(test(:,2), test(:,3),norm_data(:,1:2));

l= findobj(gcf,'tag','legend'); 
for n = 1:length(h)
   set(h(n), 'Marker', markers{unique_dates == S(n,1)})
   set(h(n), 'MarkerFaceColor', colors(unique_classes == S(n,2), :))
      set(h(n), 'MarkerEdgeColor', colors(unique_classes == S(n,2), :))       
      set(h(n), 'MarkerSize', 9);

LegendString{n} = [run_opts.date{S(n,1)} ' , ' classes{S(n,2)}];
end
set(l, 'String', LegendString)
title('PC 2 and 3')

xlabel('PC 2')
ylabel('PC 3')
xaxis = get(gca, 'xlim');
set(gca, 'xlim', [xaxis(1), xaxis(2)+2])
set(l,'location','southeast');

% figure;
% 
%      gscatter(1:size(idx,1)',idx,norm_data(:,1))
%      
%     changes =  find(diff(norm_data(:,1))~=0);
%     hold on
%     for i = 1:length(changes)
%     plot([changes(i), changes(i)], get(gca, 'ylim'), 'k--')
%     end
%     
% set(gca, 'ylim', [0 length(unique(norm_data(:,1)))+1])
% set(gca, 'xlim', [1 size(idx,1)])
% l= findobj(gcf,'tag','legend'); 
% set(l,'location','northwest');
% 
% 
% figure;
% set(gcf, 'Position', [10 10 1120 840])
% data_rows = boolean(floor(sum(~isnan(data_large),2)/size(data_large,2)));
% data_large(data_large(:,2) > 2, 2) = 3;
% 
% gscatter(data_large(data_rows,3), data_large(data_rows,4),data_large(data_rows,2))
% xlabel('RF Size (relative to parasols')
% ylabel('zc (relative to parasols')
% set(l,'location','best');

% % % set(l, 'String', {classes{unique(norm_data(:,1))}})
% data_large(data_large(:,2) == 3, 2) = 2;
% data_large(data_large(:,2) == 5, 2) = 3;
% data_large(data_large(:,2) == 6, 2) = 4;
% data_to_save = [data_large(:,1:2) test];
% csvwrite(['/Volumes/Lab/Users/crhoades/Colleen/matlab/private/colleen/New Cell Types/Cell Properties/', 'othertypestest'], data_to_save);

