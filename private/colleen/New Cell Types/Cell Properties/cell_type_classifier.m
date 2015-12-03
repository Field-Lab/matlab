% clustering cell types
clear
date = {'2005-01-21-3';'2007-01-23-5';'2007-05-01-2';'2007-08-21-4';'2009-04-13-1';'2008-04-08-2'};
concate = {'data000';'data001';'data000';'data000';'data005'; 'data003'};
    data= [];
episolon = 0.0001;
for t = 1:1%size(date,1)
    filename = ['/Volumes/Lab/Users/crhoades/Cell Properties/', date{t}, '/', concate{t}];
    classes= {'OFF Parasol', 'ON Parasol', 'ON Midget', 'OFF Midget', 'ON Large 1', 'ON Large 2', 'ON Large 3', 'ON Large 4', 'OFF Large 1', 'OFF Large 2', 'OFF Large 3', 'OFF Large 4'};
    classes= {'OFF Parasol', 'ON Parasol',  'ON Large 1', 'ON Large 2', 'ON Large 3', 'ON Large 4', 'OFF Large 1', 'OFF Large 2', 'OFF Large 3', 'OFF Large 4'};
    classes = {'ON Parasol', 'ON large 1', 'ON large 2', 'OFF Parasol', 'OFF large 1', 'OFF large 2'};
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
            new_data = [c*ones(size(cell2mat(struct2cell(output.output.parameters)'),1),1), cell2mat(struct2cell(output.output.parameters)')];
            new_data = [new_data(:,1:4), new_data(:,6:end)]; % remove t_t
            if c == 1
                ref_on_parasol = new_data;
            elseif c == 4
                ref_off_parasol = new_data;
            end
            
            if c == 2 || c == 3
                real_ind = boolean(floor(sum(~isnan(ref_on_parasol),2)/size(ref_on_parasol,2)));
                new_data(:,2:end) = new_data(:,2:end)./repmat(mean(ref_on_parasol(real_ind, 2:end)), size(new_data,1),1);
                mean_amp = mean(ref_on_parasol(real_ind, end));
                if mean_amp < 0 
                    new_data(:,end) = new_data(:,end)*-1;
                end
                
                data = [data; new_data];

            elseif c == 5 ||c ==6
                real_ind = boolean(floor(sum(~isnan(ref_off_parasol),2)/size(ref_off_parasol,2)));
                new_data(:,2:end) = new_data(:,2:end)./repmat(mean(ref_off_parasol(real_ind, 2:end)), size(new_data,1),1);
                        data = [data; new_data];
                 mean_amp = mean(ref_off_parasol(real_ind, end));
                if mean_amp < 0 
                    new_data(:,end) = new_data(:,end)*-1;
                end
            end
            
            %         new_data_unit = new_data/std(new_data);
%             data = [data; new_data];
    %         colors = [colors ;repmat(col(counter, :), size(cell2mat(struct2cell(output.output.parameters)'),1),1)];

            %         colors = [colors ;repmat(col(counter, :), size(cell2mat(struct2cell(output.output.parameters)'),1),1)];
            counter = counter+1;
        end

    end
end

data_rows = boolean(floor(sum(~isnan(data),2)/size(data,2)));
norm_data = data(:,2:end)./(repmat(std(data(data_rows,2:end)), size(data,1),1)+episolon);
norm_data = [data(:,1), norm_data];


[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(norm_data(:,2:end), 'NumComponents',2);

test = SCORE;
% idx = kmeans(norm_data(:,2:end),length(unique(norm_data(:,1))));
idx = kmeans(test,length(unique(norm_data(:,1))));
% idx = kmeans(test,4);

figure;
% subplot(3,1,1)
gscatter(test(:,1), test(:,2),norm_data(:,1))

l= findobj(gcf,'tag','legend'); 
set(l,'location','best');
set(l, 'String', {classes{unique(norm_data(:,1))}})
% subplot(3,1,2)
% gscatter(test(:,1), test(:,3),norm_data(:,1))
% subplot(3,1,3)
% gscatter(test(:,2), test(:,3),norm_data(:,1))

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