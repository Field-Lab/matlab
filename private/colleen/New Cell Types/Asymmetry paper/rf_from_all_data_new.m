% load database
clear
close all
monkey_dates = cell(0);
cell_types = [1 2 3 4];
%% read spreadsheet
[~,~,RAW]=xlsread('/Volumes/Lab/Users/crhoades/Database spreadsheet.xlsx', 'Animal');

%% Format the dates
for i = 1:size(RAW,1)
    if ~isa(RAW{i,1}, 'char') && ~isnan(RAW{i,1}) % not a date with -A/B/C or NAN
        RAW{i,1} = datestr(RAW{i,1}+datenum('30-Dec-1899'), 29);
    end
    if strcmpi(RAW{i,2}, 'Monkey')
        monkey_dates = [monkey_dates; {RAW{i,1}}];
    end
end


% Populate RAW_piece which has all the pieces, array size, ecc, and angle
[~,~,RAW_piece]=xlsread('/Volumes/Lab/Users/crhoades/Database spreadsheet.xlsx', 'Piece');
for i = 1:size(RAW_piece,1)
    array_info = RAW_piece{i,3};
    if ~isempty(strfind(num2str(array_info), '120'))
        RAW_piece{i,3} = '120';
    elseif ~isempty(strfind(num2str(array_info), '1501'))
        RAW_piece{i,3} = '30';
    elseif    ~isempty(strfind(num2str(array_info), '1504'))
        RAW_piece{i,3} = '30';
    else
        RAW_piece{i,3} = '60';
    end
    
    
    if isnumeric(RAW_piece{i,9})
        RAW_piece{i,9} = RAW_piece{i,9}*4*pi; % clock angle in radians
    else
        RAW_piece{i,9} = nan;
    end
    if ~isnumeric(RAW_piece{i,8})
        
        RAW_piece{i,8} = nan;
    end
    
    
    if ~isnan(RAW_piece{i,9})  && ~isnan(RAW_piece{i,8})
        [X,Y] = pol2cart(RAW_piece{i,9}+pi/2, RAW_piece{i,8}); % referenc epoint of all_pieces{i,3} is 12:00, change to be 3:00 by adding pi/2
        %         test = [test;all_pieces{i,3}];
        if RAW_piece{i,9} > (pi)
            E = sqrt((X*0.61)^2+Y^2); %nasal % formula in the 2002 Chichilnisky paper is wrong
            
        else
            
            E = sqrt(X^2+Y^2); %temporal
            
        end
        % replace ecc value in RAW_piece with temporal equivalent ecc
        RAW_piece{i,8} = E; % converted ecc
    else
        RAW_piece{i,8} = nan; %visual angle
        
    end
    
    
    
    
    if ~isa(RAW_piece{i,1}, 'char') && ~isnan(RAW_piece{i,1}) % not a date with -A/B/C or NAN
        RAW_piece{i,1} = datestr(RAW_piece{i,1}+datenum('30-Dec-1899'), 29);
    end
end


% all_pieces = RAW_piece(:,1);
% all_pieces(:,2) = RAW_piece(:,8); % ecc distance
% all_pieces(:,3) = RAW_piece(:,9); % angle
% all_pieces(:,5) = RAW_piece(:,3); %array id
%
% test = [];
% % all_pieces has the corrected eccentriicty measure in column 4
% for i = 2:size(all_pieces,1)
%
%
% end

%% eliminate pieces that aren't primates and put into variable monkey_pieces
results = zeros(size(monkey_dates,1), size(RAW_piece,1));
for i = 1:size(monkey_dates,1)
    for j = 1:size(RAW_piece)
        results(i,j) = ~isempty(strfind(RAW_piece{j,1}, monkey_dates{i,1}(1:10)));
    end
end

% monkey_pieces = cell(0);
to_include = sum(results);
to_delete = [];
for i =1:size(RAW_piece,1)
    if ~to_include(i) > 0
        %         if isempty(strfind(RAW_piece{i,1}, 'deleted'))
        %         to_add = strfind(RAW_piece{i,1}, ' ');
        %             if ~isempty(to_add)
        %                 monkey_pieces = [monkey_pieces; all_pieces{i,1}(1:to_add-1)];
        %     else
        to_delete = [to_delete; i];
        
        %                 RAW_piece{i,:} = cell[];
    end
    
    %     end
    
    %     end
end

% I think it removes only the first row
RAW_piece(to_delete,:) = [];

%% delete pieces that don't have analysis
to_delete = [];
monkey_pieces_with_ana = cell(0);
for i = 1:size(RAW_piece,1)
    path = ['/Volumes/Analysis/' RAW_piece{i,1}];
    e = exist(path ,'dir');
    if e ~= 7
        to_delete = [to_delete; i];
        %         monkey_pieces_with_ana = [monkey_pieces_with_ana ; monkey_pieces{i,1}];
    end
    
end
RAW_piece(to_delete,:) = [];


%% Find white noise runs for each piece
[~,~,RAW_data]=xlsread('/Volumes/Lab/Users/crhoades/Database spreadsheet.xlsx', 'Dataruns');

match = zeros(size(RAW_piece,1), size(RAW_data,1));
for i = 1:size(RAW_piece,1)
    for j = 1:size(RAW_data,1)
        match(i,j) = strcmp(RAW_data{j,1}, RAW_piece{i,1});
    end
    
end

%only search pieces that have analysis
for i = 1:size(RAW_piece,1)
    start_index = find(match(i,:) == 1);
    RAW_piece{i,12} = start_index;
    j = start_index + 1;
    while j < size(RAW_data,1)
        
        if ~isnan(RAW_data{j,1}(1))
            end_index = j-1;
            RAW_piece{i,13} = end_index;
            break;
        end
        
        
        j = j+1;
    end
    
end



%% get all stixel widths

% parameters is somethingx4 for each piece where the datarun, interval,
% stixel_size, and type of display are indexed for each run
for i = 1:size(RAW_piece,1)
    if ~isempty(RAW_piece{i,12}) && ~isempty(RAW_piece{i,13})
        for j = 1:RAW_piece{i,13} - RAW_piece{i,12}+1;
            if strcmp(RAW_data{RAW_piece{i,12}+j-1, 13}, 'binary')
                parameters{i,1}{j,1} = RAW_data{RAW_piece{i,12}+j-1, 4}; %datarun
                parameters{i,1}{j,2} = RAW_data{RAW_piece{i,12}+j-1, 15}; % interval
                
                parameters{i,1}{j,3} = RAW_data{RAW_piece{i,12}+j-1, 19}; % stixel_size
                parameters{i,1}{j,4} = RAW_data{RAW_piece{i,12}+j-1, 8}; % display
                parameters{i,1}{j,5} = RAW_piece{i,8}; % equi ecc value
                parameters{i,1}{j,6} = RAW_data{RAW_piece{i,12}+j-1, 6}; % ndf data
                parameters{i,1}{j,7} = RAW_data{RAW_piece{i,12}+j-1, 10}; % objective
                parameters{i,1}{j,8} = RAW_piece{i,3}; % array size
                parameters{i,1}{j,9} = RAW_data{RAW_piece{i,12}+j-1, 14}; % RGB or BW
                parameters{i,1}{j,10} = RAW_data{RAW_piece{i,12}+j-1, 17}; % seed
                parameters{i,1}{j,11} = RAW_data{RAW_piece{i,12}+j-1, 31}; % jitter
                
            else
                parameters{i,1}{j,1} = nan;
                parameters{i,1}{j,2} = nan;
                parameters{i,1}{j,3} = nan;
                parameters{i,1}{j,4} = nan;
                parameters{i,1}{j,5} = nan;
                parameters{i,1}{j,6} = nan;
                parameters{i,1}{j,7} = nan;
                parameters{i,1}{j,8} = nan;
                parameters{i,1}{j,9} = nan;
                parameters{i,1}{j,10} = nan;
                parameters{i,1}{j,11} = nan;
                
            end
            
        end
    end
end

% find out if ON and OFF midgets are in analysis
good_data_sets = cell(0);
counter = 1;
pieces_to_include = cell(0);
should_run = cell(0);

for i = 1:size(parameters,1)
    match2 = zeros(size(parameters,1), size(RAW_piece,1));
    
    clear datarun
    clear param_choices
    
    test = [];
    for j = 1:size(parameters{i,1},1)
        test(j) = ~isnan(parameters{i,1}{j,1}(1));
    end
    
    if ~isempty(test) && sum(test) > 0
        %         if sum(test)>1
        possibilities = find(test == 1);
        for t = possibilities
            
            
            param_choices(t,1) = parameters{i,1}{t,3}; % stixel size
            if strcmp(parameters{i,1}{t,4}, 'CRT')
                
                param_choices(t,2) = 1; % CRT
            elseif strcmp(parameters{i,1}{t,4}, 'OLED')
                param_choices(t,2) = 2; % OLED
            else
                param_choices(t,2) = 0; % OLED
                
            end
            param_choices(t,3) = parameters{i,1}{t,5}; %ecc
            
            if ischar(parameters{i,1}{t,6})
                param_choices(t,4) = 0;
            else
                param_choices(t,4) = parameters{i,1}{t,6}; % NDF
            end
            
            if ischar(parameters{i,1}{t,7})
                param_choices(t,5) = 0;
            else
                param_choices(t,5) = parameters{i,1}{t,7}; % objective
            end
            
            
            
            param_choices(t,6) = str2double(parameters{i,1}{t,8}); % array size
            if strcmp(parameters{i,1}{t,9}, 'BW')
                param_choices(t,7) = 1; % BW
                
            else
                param_choices(t,7) = 3; % RGB
            end
            
            if ischar(parameters{i,1}{t,10})
                param_choices(t,8) = nan;
            else
                
                param_choices(t,8) = parameters{i,1}{t,10}; % Seed
            end
            
            param_choices(t,9) = parameters{i,1}{t,2}; % Interval
            
            
            
            if strcmp(parameters{i,1}{t,11}, 'T')
                param_choices(t,10) = 1;
            else
                
                param_choices(t,10) = 0; % Seed
            end
            
        end
        
        for p = 1:size(param_choices,1)
            
            %             if param_choices(p,9) <4 % interval
            %                 param_choices(p,:) = inf;
            %             end
            
            if param_choices(p,2) == 2 % OLED
                param_choices(p,:) = inf;
                param_choices(p,:) = inf;
            end
            
            if param_choices(p,3) < 12 % ecc < 10mm
                param_choices(param_choices(:,1) <= 1,:) = inf; %implement stixel_choices based on ecc?
                param_choices(param_choices(:,1) >= 10,:) = inf;
            else % ecc >= 8
                param_choices(param_choices(:,1) <= 2,:) = inf; %implement stixel_choices based on ecc?
                param_choices(param_choices(:,1) >= 12,:) = inf;
            end
            
            if param_choices(p,4) > 0.6 % NDF larger than 0.6
                param_choices(p,:) = inf;
            end
            
            if param_choices(p,5) < 6.5 &&  param_choices(p,5) ~= 0 % objective only 6.5
                param_choices(p,:) = inf;
            end
            %
            %              if param_choices(p,5) == 2 % obejective only 6.5 and 4
            %                 param_choices(p,:) = inf;
            %             end
            
            if param_choices(p,5) == 10 % obejective not 10
                param_choices(p,:) = inf;
            end
            
            if param_choices(p,6) > 60 || param_choices(p,6) < 60 % array not 60um
                param_choices(p,:) = inf;
            end
            if param_choices(p,10) == 1 % jitter
                param_choices(p,:) = inf;
            end
        end
        
        %             opt_stixel_choices = find(~isinf(stixel_choices));
        %             [stixel_size, smallest_available] = min(stixel_choices);
        %             interval = smallest_available;
        for o = 1:size(param_choices,1)
            try
                if ~isinf(param_choices(o, 1)) && ~isinf(param_choices(o, 2))
                    if exist(['/Volumes/Analysis/', RAW_piece{i,1}, '/', parameters{i,1}{o,1}, '/',parameters{i,1}{o,1},'.params'], 'file') ~= 0
                        
                        datarun = load_data([RAW_piece{i,1}, '/', parameters{i,1}{o,1}]);
                        datarun = load_params(datarun);
                        
                        
                        %                         if size(datarun.cell_types{1}.cell_ids,2) > 0 && size(datarun.cell_types{2}.cell_ids,2) > 0 && size(datarun.cell_types{3}.cell_ids,2) > 0 && size(datarun.cell_types{4}.cell_ids,2) > 0
                        for z = 1:length(cell_types) % has 10 of all targeted cell types
                            if size(datarun.cell_types{cell_types(z)}.cell_ids,2) > 10
                                % do nothing
                            else
                                
                                param_choices(o,:) = inf;
                                
                            end
                        end
                    else
                        param_choices(o,:) = inf;
                        
                        if strcmp(parameters{i,1}{o,9}, 'RGB')
                            type = 'RGB';
                        else
                            type = 'BW';
                        end
                        stimulus = [type, '-', num2str(parameters{i,1}{o,3}), '-',num2str(parameters{i,1}{o,2}), '-0.48-',num2str(parameters{i,1}{o,10}), '.xml'];
                        if parameters{i,1}{o,2}<=6 % interval <=6
                            should_run = [should_run; {RAW_piece{i,1}, parameters{i,1}{o,1}, stimulus}];
                        end
                    end
                    
                end
            catch
                disp([RAW_piece{i,1}, '/', parameters{i,1}{o,1}])
                
            end
            
        end
        if sum(param_choices(param_choices(:,1)<inf,1))>0 % there is an option
            %                 continue;
            %             else
            poss = find(param_choices(:,1)<inf);
            [stixel_size, ind] = min(param_choices(poss,1));
            ecc = param_choices(poss(ind),3);
            
            interval = poss(ind);
        else
            skip(i) = 1;
            continue;
        end
        
        
        k = interval;
        
        try
            datarun = load_data([RAW_piece{i,1}, '/', parameters{i,1}{k,1}]);
            datarun = load_params(datarun);
            datarun = load_globals(datarun);
            datarun = load_ei(datarun, 'all');
            
            
            %                 datarun = load_sta(datarun, 'load_sta', datarun.cell_ids(1), 'verbose', 0);
            
            %         if ~strcmp(cell2mat([RAW_piece(i,1), parameters{i,1}{k,1}]), '2008-04-30-2data007') && ~strcmp(cell2mat([RAW_piece(i,1), parameters{i,1}{k,1}]), '2012-04-13-5data001') && ~strcmp(cell2mat([RAW_piece(i,1), parameters{i,1}{k,1}]), '2008-07-07-3data004')% this one is messed up RF sizes
            RAW_piece(i,14:14+size(parameters{i,1}(k,:),2)-1) = parameters{i,1}(k,:); % add the datarun number picked
            pieces_to_include = [pieces_to_include ; RAW_piece(i,:)];
            
            for m = 1:size(datarun.vision.sta_fits,1)
                x_mean(m) = datarun.vision.sta_fits{m}.mean(1);
                y_mean(m) = datarun.vision.sta_fits{m}.mean(2);
            end
            xstart = prctile(x_mean,10);
            xend = prctile(x_mean,90);
            rf_area = xend-xstart; %rf_area in stixel units,  should be approx 2mm
            %
            
            %             center_spacing
            
            %             if isnan(param_choices(k,5))
            %                 param_choices(k,5) = 6.5;
            %             end
            
            
            for r = cell_types
                cell_index = get_cell_indices(datarun, datarun.cell_types{r}.cell_ids);
                cell_type_mean1 = [];
                for c = 1:length(cell_index)
                    cell_type_mean1(c,1) = datarun.vision.sta_fits{cell_index(c)}.mean(1);
                    cell_type_mean1(c,2) = datarun.vision.sta_fits{cell_index(c)}.mean(2);
                    
                    
                    
                    rfs{counter}{r}(c) = sqrt(datarun.vision.sta_fits{cell_index(c)}.sd(1)*datarun.vision.sta_fits{cell_index(c)}.sd(2))*2*RAW_piece{i,16}*5; %2 radii/diameter*pixels/stixel*5um/pixel
                    tc_r{counter}{r}(c,:) =  datarun.vision.timecourses(cell_index(c)).r;
                    tc_g{counter}{r}(c,:) =  datarun.vision.timecourses(cell_index(c)).g;
                    tc_b{counter}{r}(c,:) =  datarun.vision.timecourses(cell_index(c)).b;
                    
                    
                    %                     rfs{counter}{r}(c) = sqrt(datarun.vision.sta_fits{cell_index(c)}.sd(1)*datarun.vision.sta_fits{cell_index(c)}.sd(2))*2*2/rf_area*1000; %rf_area/2mm*1000um/mm
                    
                    %                                                         rfs{counter}{r}(c) = sqrt(datarun.vision.sta_fits{cell_index(c)}.sd(1)*datarun.vision.sta_fits{cell_index(c)}.sd(2))*2*stixel_size*5*6.5/scale;
                    
                end
                %                 D = pdist(cell_type_mean1);
                spacing(counter,r) = ei_center(datarun,cell_index);
                
                [~,d] = knnsearch(cell_type_mean1, cell_type_mean1, 'K', 5, 'IncludeTies', true);
                me = zeros(size(d,1),1);
                for ii = 1:size(d,1)
                    me(ii) = mean(d{ii}(2:end));
                end
                inter_cell_distance(counter, r) =median(me)*RAW_piece{i,16}*5; %units of um
                tc{counter}{r}(1,:) = median(tc_r{counter}{r});
                tc{counter}{r}(2,:) = median(tc_g{counter}{r});
                tc{counter}{r}(3,:) = median(tc_b{counter}{r});
                
            end
            counter= counter +1;
            %         end
            
        catch
                disp([RAW_piece{i,1}, '/', parameters{i,1}{o,1}])
        end
    end
end


% Write the "should write" to csv file
% xlswrite('to_run', should_run);

%% Process TCs

for i = cell_types
    for j = 1:size(tc,2)
        [a, max_ind] = max(tc{j}{i}(2,:)); % max of green channel
        [b, min_ind] = min(tc{j}{i}(2,:)); % min of green channel
        biphasic(j,i) = max(abs(a),abs(b))/ min(abs(a),abs(b));
        
        
        if max_ind > min_ind
            peak_ind = max_ind;
        else
            peak_ind = min_ind;
        end
        
        [t_ind,t0] = crossing(tc{j}{i}(2,:));
        
        less_than = find(t_ind < peak_ind == 1);
        if isempty(less_than)
            t_zc = nan;
        else
            
            ind = less_than(end);
            % [~,ind] = min(abs(t0 - run_opts.num_frames/2)); % find the frame that closest to the middle of num_frames
            t0 = t0(ind);
            refresh =pieces_to_include{j,15}/120*1000;
            t_zc(j,i) = abs((30-t0)*refresh); % 30 frames in STA
        end
        
        
        
    end
end


% Parasol biphasic index
on_parasol_increase = biphasic(biphasic(:,1)<7,1)./biphasic(biphasic(:,1)<7,2);
pec_inc_parasol = median(on_parasol_increase);

figure;
plot(biphasic(biphasic(:,1)<7,1),biphasic(biphasic(:,1)<7,2), 'ok', 'MarkerFaceColor', 'k')
hold on
xlim = get(gca, 'xlim');
ylim = get(gca,'ylim');
plot([0 max(xlim(2), ylim(2))], [0 max(xlim(2), ylim(2))], '--k')
xlabel('ON parasol Biphasic Index')
ylabel('OFF parasol Biphasic Index')
title(['ON parasols have a biphasic index that is on average ', num2str(round(pec_inc_parasol*1000)/10-100), '% higher'])

% midget biphasic index
on_midget_increase = biphasic(:,4)./biphasic(:,3);
pec_inc_midget = median(on_midget_increase);

figure;
plot(biphasic(:,3),biphasic(:,4), 'ok', 'MarkerFaceColor', 'k')
hold on
xlim = get(gca, 'xlim');
ylim = get(gca,'ylim');
plot([0 max(xlim(2), ylim(2))], [0 max(xlim(2), ylim(2))], '--k')
xlabel('ON midget Biphasic Index')
ylabel('OFF midget Biphasic Index')
title(['OFF midgets have a biphasic index that is on average ', num2str(round(pec_inc_midget*1000)/10-100), '% higher'])




% Parasol zero crossing
on_parasol_increase = t_zc(t_zc(:,1)<150,2)./t_zc(t_zc(:,1)<150,1);
pec_inc_parasol = median(on_parasol_increase);

figure;
plot(t_zc(t_zc(:,1)<150,1),t_zc(t_zc(:,1)<150,2), 'ok', 'MarkerFaceColor', 'k')
hold on
xlim = get(gca, 'xlim');
ylim = get(gca,'ylim');
plot([0 max(xlim(2), ylim(2))], [0 max(xlim(2), ylim(2))], '--k')
xlabel('ON parasol Zero Crossing (ms)')
ylabel('OFF parasol Zero Crossing (ms)')
title(['ON parasols are faster than OFF parasols on average by ', num2str(round(pec_inc_parasol*1000)/10-100), '%'])

% midget zero crossing
on_midget_increase = t_zc(t_zc(:,3)<200,4)./t_zc(t_zc(:,3)<200,3);
pec_inc_midget = median(on_midget_increase);

figure;
plot(t_zc(t_zc(:,3)<200,3),t_zc(t_zc(:,3)<200,4), 'ok', 'MarkerFaceColor', 'k')
hold on
xlim = get(gca, 'xlim');
ylim = get(gca,'ylim');
plot([0 max(xlim(2), ylim(2))], [0 max(xlim(2), ylim(2))], '--k')
xlabel('ON midget Zero Crossing (ms)')
ylabel('OFF midget Zero Crossing (ms)')
title(['ON midgets are faster than OFF midgets on average by ', num2str(round(pec_inc_midget*1000)/10-100), '%'])



%% Process RFs
for i = cell_types
    for j= 1:size(rfs,2)
        if size(rfs{1,j}{i},2) > 10
            [B, IDX, OUTLIERS] = deleteoutliers(rfs{1,j}{i}, 0.05, 1);
            
            rf_cell_type(j,i) = median(B(~isnan(B)));
        else
            rf_cell_type(j,i) = nan;
        end
        
    end
end

% rf_cell_type(43,:) = nan;
rf_cell_type_include = rf_cell_type(~any(isnan(rf_cell_type),2),:);
on_parasol_increase = rf_cell_type_include(:,1)./rf_cell_type_include(:,2);
pec_inc_parasol = median(on_parasol_increase(~isnan(on_parasol_increase)))
on_midget_increase =rf_cell_type_include(:,3)./rf_cell_type_include(:,4);
pec_inc_midget = median(on_midget_increase(~isnan(on_midget_increase)))

figure; plot(rf_cell_type(:,1), rf_cell_type(:,2), 'ok', 'MarkerFaceColor', 'k')
hold on
xlim = get(gca, 'xlim');
ylim = get(gca,'ylim');
plot([0 max(xlim(2), ylim(2))], [0 max(xlim(2), ylim(2))], '--k')
xlabel('ON parasol RF diameter (\mum)')
ylabel('OFF parasol RF diameter (\mum)')
axis equal
axis square
title(['ON parasols are larger than OFF parasols on average by ', num2str(round(pec_inc_parasol*1000)/10-100), '%'])


figure; plot(rf_cell_type(:,3), rf_cell_type(:,4), 'ok', 'MarkerFaceColor', 'k')
hold on
xlim = get(gca, 'xlim');
ylim = get(gca,'ylim');
plot([0 max(xlim(2), ylim(2))], [0 max(xlim(2), ylim(2))], '--k')
xlabel('ON midget RF diameter (\mum)')
ylabel('OFF midget RF diameter (\mum)')
title(['ON midgets are larger than OFF midgets on average by ', num2str(round(pec_inc_midget*1000)/10-100), '%'])
axis equal

% for i = 1:size(good_data_sets(:,6),1)
%     if isnumeric(good_data_sets{i,6}) == 0
%         ecc(i) = 0;
%     else
%         ecc(i) = good_data_sets{i,6};
%     end
%
% end

% ecc_include = ecc(~any(isnan(rf_cell_type),2)');
% figure; gscatter(rf_cell_type_include(~isnan(ecc_include),1), rf_cell_type_include(~isnan(ecc_include),2), ecc_include(~isnan(ecc_include)))
% hold on
% xlim = get(gca, 'xlim');
% ylim = get(gca,'ylim');
% plot([0 max(xlim(2), ylim(2))], [0 max(xlim(2), ylim(2))], '--k')
% xlabel('ON parasol')
% ylabel('OFF parasol')

figure;
% plot(ecc_include(~isnan(ecc_include)), rf_cell_type_include(~isnan(ecc_include),1), 'ko', 'MarkerFaceColor', 'k')
hold on
plot(cell2mat(pieces_to_include(rf_cell_type(:,1)>20, 18)), rf_cell_type(rf_cell_type(:,1)>20,1), 'ko', 'MarkerFaceColor', 'w')
plot(cell2mat(pieces_to_include(rf_cell_type(:,1)>20, 18)), rf_cell_type(rf_cell_type(:,1)>20,2), 'ko', 'MarkerFaceColor', 'k')

xlabel('Temporal equivalent eccentricity (mm)')
ylabel('RF diameter (\mum)')
title('Parasols')
legend('ON parasol', 'OFF parasol');

% for i =1:size(good_data_sets,1)
%     if rf_cell_type(i,3)<80 && ecc(i) < 8
%         rf_cell_type(i,3:4) = nan;
%     end
% end

figure;
% plot(ecc_include(~isnan(ecc_include)), rf_cell_type_include(~isnan(ecc_include),1), 'ko', 'MarkerFaceColor', 'k')
hold on
plot(cell2mat(pieces_to_include(:, 18)), rf_cell_type(:,3), 'ko', 'MarkerFaceColor', 'w')
plot(cell2mat(pieces_to_include(:, 18)), rf_cell_type(:,4), 'ko', 'MarkerFaceColor', 'k')

xlabel('Temporal equivalent eccentricity (mm)')
ylabel('RF diameter (\mum)')
title('Midgets')
legend('ON midget', 'OFF midget');


% b1 = [ones(length(ecc_include(~isnan(ecc_include))'),1) ecc_include(~isnan(ecc_include))']\rf_cell_type_include(~isnan(ecc_include),2);
% yCalc1 = [ones(length(ecc_include(~isnan(ecc_include))'),1) ecc_include(~isnan(ecc_include))']*b1;
% hold on
% plot(ecc_include(~isnan(ecc_include)), yCalc1, 'o')

%
% figure; plot(ecc_include(~isnan(ecc_include)), rf_cell_type_include(~isnan(ecc_include),3), 'ko', 'MarkerFaceColor', 'k')
% hold on
% plot(ecc_include(~isnan(ecc_include)), rf_cell_type_include(~isnan(ecc_include),4), 'ro', 'MarkerFaceColor', 'r')
% title('Midget versus ecc')



% figure; gscatter(rf_cell_type_include(~isnan(ecc_include),3), rf_cell_type_include(~isnan(ecc_include),4), ecc_include(~isnan(ecc_include)))
% hold on
% xlim = get(gca, 'xlim');
% ylim = get(gca,'ylim');
% plot([0 max(xlim(2), ylim(2))], [0 max(xlim(2), ylim(2))], '--k')
% xlabel('ON midget')
% ylabel('OFF midget')


% %% Inter cell distance

figure; plot(cell2mat(pieces_to_include(:, 18)), inter_cell_distance(:,1), 'ko', 'markerfacecolor', 'k') % on parasol intracell distance
xlabel('Temporal equivalent eccentricity')
ylabel('STA Center to Center Spacing (ON parasols in \mum)')
title('Variability in STA center-to-center spacing')

figure; plot(cell2mat(pieces_to_include(:, 18)), inter_cell_distance(:,2), 'ko', 'markerfacecolor', 'k') % on parasol intracell distance
xlabel('Temporal equivalent eccentricity')
ylabel('STA Center to Center Spacing (OFF parasols in \mum)')
title('Variability in STA center-to-center spacing')


figure; plot(inter_cell_distance(:,1), rf_cell_type(:,1), 'ko', 'markerfacecolor', 'k') % on parasol intracell distance
hold on 
xaxis = get(gca, 'xlim');
yaxis = get(gca, 'ylim');
plot([min(xaxis(1), yaxis(1)), max(xaxis(2), yaxis(2))], [ min(xaxis(1), yaxis(1)) , max(xaxis(2), yaxis(2))], '--k')
ylabel('RF Diameter (\mum)')
xlabel('STA Center to Center Spacing (ON parasols in \mum)')
title('Mosaic Spacing Correlated with RF Size')


figure; plot(spacing(:,1), rf_cell_type(:,1), 'ko', 'markerfacecolor', 'k') % on parasol intracell distance
hold on 
xaxis = get(gca, 'xlim');
yaxis = get(gca, 'ylim');
plot([min(xaxis(1), yaxis(1)), max(xaxis(2), yaxis(2))], [ min(xaxis(1), yaxis(1)) , max(xaxis(2), yaxis(2))], '--k')
ylabel('RF Diameter (\mum)')
xlabel('EI Center to Center Spacing (ON parasols in \mum)')
title('Mosaic Spacing Correlated with RF Size')


figure; plot(spacing(:,1), inter_cell_distance(:,1),'ko', 'markerfacecolor', 'k')
hold on 
xaxis = get(gca, 'xlim');
yaxis = get(gca, 'ylim');
plot([min(xaxis(1), yaxis(1)), max(xaxis(2), yaxis(2))], [ min(xaxis(1), yaxis(1)) , max(xaxis(2), yaxis(2))], '--k')

title('ON parasol')
xlabel('EI center to center spacing (\mum)')
ylabel('STA center to center spacing (\mum)')

scale_factor_on = inter_cell_distance(:,1)./spacing(:,1);
scale_factor_off = inter_cell_distance(:,2)./spacing(:,2);

figure;
hold on
plot(cell2mat(pieces_to_include(rf_cell_type(:,1)>20, 18)), rf_cell_type(rf_cell_type(:,1)>20,1)./scale_factor_on, 'ko', 'MarkerFaceColor', 'w')
plot(cell2mat(pieces_to_include(rf_cell_type(:,1)>20, 18)), rf_cell_type(rf_cell_type(:,1)>20,2)./scale_factor_off, 'ko', 'MarkerFaceColor', 'k')
xlabel('Temporal equivalent eccentricity (mm)')
ylabel('RF diameter (\mum)')
title('Parasols')
legend('ON parasol', 'OFF parasol');

figure; plot(spacing(:,2), inter_cell_distance(:,2), 'ko', 'markerfacecolor', 'k')
hold on 
xaxis = get(gca, 'xlim');
yaxis = get(gca, 'ylim');
plot([min(xaxis(1), yaxis(1)), max(xaxis(2), yaxis(2))], [ min(xaxis(1), yaxis(1)) , max(xaxis(2), yaxis(2))], '--k')
title('OFF parasol')
xlabel('EI center to center spacing (\mum)')
ylabel('STA center to center spacing (\mum)')


figure; plot(cell2mat(pieces_to_include(:, 18)), spacing(:,1), 'ko', 'markerfacecolor', 'k') % on parasol intracell distance
xlabel('Temporal equivalent eccentricity')
ylabel('EI Center to Center Spacing (ON parasols in \mum)')
title('Variability in EI center-to-center spacing')

figure; plot(cell2mat(pieces_to_include(:, 18)), spacing(:,2), 'ko', 'markerfacecolor', 'k') % on parasol intracell distance
xlabel('Temporal equivalent eccentricity')
ylabel('EI Center to Center Spacing (OFF parasols in \mum)')
title('Variability in EI center-to-center spacing')



[~,git_hash_string] = system('git rev-parse HEAD');
fprintf('Git Hash: %s \n', git_hash_string);

