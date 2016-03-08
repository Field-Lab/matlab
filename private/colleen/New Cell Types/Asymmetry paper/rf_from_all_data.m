% load database
clear
monkey_dates = cell(0);

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
    
    
    if ~isa(RAW_piece{i,1}, 'char') && ~isnan(RAW_piece{i,1}) % not a date with -A/B/C or NAN
        RAW_piece{i,1} = datestr(RAW_piece{i,1}+datenum('30-Dec-1899'), 29);
    end
end


all_pieces = RAW_piece(:,1);
all_pieces(:,2) = RAW_piece(:,8); % ecc distance
all_pieces(:,3) = RAW_piece(:,9); % angle
all_pieces(:,5) = RAW_piece(:,3); %array id

test = [];
% all_pieces has the corrected eccentriicty measure in column 4
for i = 2:size(all_pieces,1)
    
    if ~isnan(all_pieces{i,3})  && ~isnan(all_pieces{i,2})
        [X,Y] = pol2cart(all_pieces{i,3}+pi/2, all_pieces{i,2}); % referenc epoint of all_pieces{i,3} is 12:00, change to be 3:00 by adding pi/2
        test = [test;all_pieces{i,3}];
        if all_pieces{i,3} > (pi)
            E = sqrt((X*0.61)^2+Y^2); %nasal % formula in the 2002 Chichilnisky paper is wrong

        else

            E = sqrt(X^2+Y^2); %temporal
            
        end
        %     all_pieces{i,4} = 0.1+4.21*E+0.038*E^2; %visual angle
        all_pieces{i,4} = E; % converted ecc
    else
        all_pieces{i,4} = nan; %visual angle
        
    end
    
end

%% eliminate pieces that aren't primates and put into variable monkey_pieces
results = zeros(size(monkey_dates,1), size(all_pieces,1));
for i = 1:size(monkey_dates,1)
    for j = 1:size(all_pieces)
        results(i,j) = ~isempty(strfind(all_pieces{j,1}, monkey_dates{i,1}(1:10)));
    end
end

monkey_pieces = cell(0);
to_include = sum(results);
for i =1:size(all_pieces,1)
    if to_include(i) > 0
        if isempty(strfind(all_pieces{i,1}, 'deleted'))
            to_add = strfind(all_pieces{i,1}, ' ');
            if ~isempty(to_add)
                monkey_pieces = [monkey_pieces; all_pieces{i,1}(1:to_add-1)];
            else
                monkey_pieces = [monkey_pieces; all_pieces{i,1}];
            end
            
        end
        
    end
end


%% find pieces with some analysis done
monkey_pieces_with_ana = cell(0);
for i = 1:size(monkey_pieces,1)
    path = ['/Volumes/Analysis/' monkey_pieces{i,1}];
    e = exist(path ,'dir');
    if e == 7
        monkey_pieces_with_ana = [monkey_pieces_with_ana ; monkey_pieces{i,1}];
    end
    
end

%% Find white noise runs for each piece
[~,~,RAW_data]=xlsread('/Volumes/Lab/Users/crhoades/Database spreadsheet.xlsx', 'Dataruns');

match = zeros(size(monkey_pieces_with_ana,1), size(RAW_data,1));
for i = 1:size(monkey_pieces_with_ana,1)
    for j = 1:size(RAW_data,1)
        match(i,j) = strcmp(RAW_data{j,1}, monkey_pieces_with_ana{i,1});
    end
    
end

%only search pieces that have analysis
for i = 1:size(monkey_pieces_with_ana,1)
    start_index = find(match(i,:) == 1);
    monkey_pieces_with_ana{i,2} = start_index;
    j = start_index + 1;
    while j < size(RAW_data,1)
        
        if ~isnan(RAW_data{j,1}(1))
            end_index = j-1;
            monkey_pieces_with_ana{i,3} = end_index;
            break;
        end
        
        
        j = j+1;
    end
    
end



%% get all stixel widths

for i = 1:size(monkey_pieces_with_ana,1)
    if ~isempty(monkey_pieces_with_ana{i,2}) && ~isempty(monkey_pieces_with_ana{i,3})
        for j = 1:monkey_pieces_with_ana{i,3} - monkey_pieces_with_ana{i,2}+1;
            if strcmp(RAW_data{monkey_pieces_with_ana{i,2}+j-1, 13}, 'binary')
                parameters{i,1}{j,1} = RAW_data{monkey_pieces_with_ana{i,2}+j-1, 4};
                parameters{i,1}{j,2} = RAW_data{monkey_pieces_with_ana{i,2}+j-1, 15}; % interval
                
                parameters{i,1}{j,3} = RAW_data{monkey_pieces_with_ana{i,2}+j-1, 19}; % stixel_size
                parameters{i,1}{j,4} = RAW_data{monkey_pieces_with_ana{i,2}+j-1, 8}; % display
                
                
            else
                parameters{i,1}{j,1} = nan;
                parameters{i,1}{j,2} = nan;
                parameters{i,1}{j,3} = nan;
                parameters{i,1}{j,4} = nan;
                
            end
            
        end
    end
end

% find out if ON and OFF midgets are in analysis
good_data_sets = cell(0);
counter = 0;
for i = 1:size(parameters,1)
    match2 = zeros(size(parameters,1), size(all_pieces,1));
    
    clear datarun
    clear stixel_choices
    clear x_mean
    clear y_mean
    test = [];
    for j = 1:size(parameters{i,1},1)
        test(j) = ~isnan(parameters{i,1}{j,1}(1));
    end
    
    if ~isempty(test) && sum(test) > 0
        if sum(test)>1
            possibilities = find(test == 1);
            for t = possibilities
                stixel_choices(t) = parameters{i,1}{t,3};
            end
            stixel_choices(stixel_choices <= 2) = inf; %implement stixel_choices based on ecc?
            stixel_choices(stixel_choices >= 12) = inf;
            
            %             opt_stixel_choices = find(~isinf(stixel_choices));
            %             [stixel_size, smallest_available] = min(stixel_choices);
            %             interval = smallest_available;
            for o = 1:length(stixel_choices)
                try
                    if ~isinf(stixel_choices(o))
                        if exist(['/Volumes/Analysis/', monkey_pieces_with_ana{i,1}, '/', parameters{i,1}{o,1}, '/',parameters{i,1}{o,1},'.params'], 'file') ~= 0
                            
                            datarun = load_data([monkey_pieces_with_ana{i,1}, '/', parameters{i,1}{o,1}]);
                            datarun = load_params(datarun);
                            
                            
                            if size(datarun.cell_types{1}.cell_ids,2) > 0 && size(datarun.cell_types{2}.cell_ids,2) > 0 && size(datarun.cell_types{3}.cell_ids,2) > 0 && size(datarun.cell_types{4}.cell_ids,2) > 0
                                % do nothing
                            else
                                
                                stixel_choices(o) = inf;
                            end
                            
                        else
                            stixel_choices(o) = inf;
                            
                        end
                        
                    end
                catch
                    disp('error loading params')
                end
                
            end
            if sum(stixel_choices)<1
                continue;
            else
                poss = find(stixel_choices<inf);
                [stixel_size, ind] = min(stixel_choices(poss));
                interval = poss(ind);
            end
            
        else
            
            interval = find(test == 1);
            stixel_size = parameters{i,1}{interval,3};
            
        end
        
        k = interval;
        %             datarun.names.rrs_params_path = ;
        
        if exist(['/Volumes/Analysis/', monkey_pieces_with_ana{i,1}, '/', parameters{i,1}{k,1}, '/',parameters{i,1}{k,1},'.params'], 'file') > 0
            %                 opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
            if ~strcmp(parameters{i,1}{k,4}, 'OLED')
                
                try
                    % Find the eccentricity
                    for co = 1:size(all_pieces,1)
                        if isempty(strfind(all_pieces{co,1}, monkey_pieces_with_ana{i,1}))
                            match2(i,co) = 0;
                        else match2(i,co) = 1;
                            break;
                        end
                        
                    end
                    
                    
                    datarun = load_data([monkey_pieces_with_ana{i,1}, '/', parameters{i,1}{k,1}]);
                    datarun = load_params(datarun);
                    %                 datarun = load_sta(datarun, 'load_sta', datarun.cell_ids(1), 'verbose', 0);
                    
                    if size(datarun.cell_types{1}.cell_ids,2) > 0 && size(datarun.cell_types{2}.cell_ids,2) > 0 && size(datarun.cell_types{3}.cell_ids,2) > 0 && size(datarun.cell_types{4}.cell_ids,2) > 0
                        if ~strcmp(cell2mat([monkey_pieces_with_ana(i,1), parameters{i,1}{k,1}]), '2008-04-30-2data007') && ~strcmp(cell2mat([monkey_pieces_with_ana(i,1), parameters{i,1}{k,1}]), '2012-04-13-5data001') && ~strcmp(cell2mat([monkey_pieces_with_ana(i,1), parameters{i,1}{k,1}]), '2008-07-07-3data004')% this one is messed up RF sizes
                            if ~strcmp(all_pieces{co,5}, '60')
                                continue;
                            end
                            b = 0;
                            monkey_pieces_with_ana{i,4} = k;
                            for ii = 1:size(RAW_data,1)
                               if ~isempty(strfind(RAW_data{ii,1},monkey_pieces_with_ana{i,1}))
                                   if ~isempty(strfind(RAW_data{ii,4}, parameters{i,1}{k,1}))
                                         ndf_data = RAW_data{ii,6}; 
                                         b = 1;

                                       end
                                       ii = ii+1;
                                   while ~strcmp(RAW_data{ii,4}, 'data000') && b == 0
                                       if ~isempty(strfind(RAW_data{ii,4}, parameters{i,1}{k,1}))
                                         ndf_data = RAW_data{ii,6}; 
                                         b = 1;

                                       end
                                       ii = ii+1;
                                   end
                                   
                               end
                                if b == 1
                                    break;
                                end
                                
                            end
                            if isnan(ndf_data) | ndf_data==0 | ndf_data==0.3| ndf_data==0.6| isempty(ndf_data) |strcmp(ndf_data, 'none')
                            else
                                
                                continue;
                            end
                            
                            good_data_sets = [good_data_sets ; [monkey_pieces_with_ana(i,1), parameters{i,1}{k,1}, monkey_pieces_with_ana(i,2:end), all_pieces{co,4}, stixel_size, all_pieces{co,5}]]; % end with objective,  ecc and stixel size and array info
                            counter = counter +1;
                            
                            for ii = 1:size(RAW_data,1)
                                if ~isempty(strfind(RAW_data{ii,1},good_data_sets{counter,1}))
                                    location = ii;
                                    break;
                                end
                            end
                            for iii = location:size(RAW_data,1)
                                if ~isempty(strfind(RAW_data{iii,4},good_data_sets{counter,2}))
                                    if strcmp(RAW_data{iii,10}, 'missing')
                                        objective(counter) = nan;
                                    else
                                        
                                        objective(counter) = RAW_data{iii,10};
                                    end
                                    
                                    break;
                                end
                                
                            end
                            
                                 if isnan(objective(counter))
                                        if strcmp( good_data_sets{counter,8}, '120')
                                            scale = 4;
                                        elseif strcmp( good_data_sets{counter,8}, '30')
                                            scale = 6.5;
                                        else
                                            scale = 6.5;
                                        end
                                 else
                                        if objective(counter) == 4
                                            objective(counter) = 100;
                                        elseif objective(counter) == 10
                                           objective(counter) = 100;

                                        end
                                        
                                        scale = objective(counter);
                                 end
                                    
                            for m = 1:size(datarun.vision.sta_fits,1)  
                                x_mean(m) = datarun.vision.sta_fits{m}.mean(1);
                                y_mean(m) = datarun.vision.sta_fits{m}.mean(2);
                            end
                            xstart = prctile(x_mean,10);
                            xend = prctile(x_mean,90);
                            rf_area = xend-xstart; %rf_area in stixel units,  should be approx 2mm

                            for r = 1:4
                                cell_index = get_cell_indices(datarun, datarun.cell_types{r}.cell_ids);
                                for c = 1:length(cell_index)
                                    
                               
                                        
                                   rfs{counter}{r}(c) = sqrt(datarun.vision.sta_fits{cell_index(c)}.sd(1)*datarun.vision.sta_fits{cell_index(c)}.sd(2))*2*2/rf_area*1000*6.5/scale; %rf_area/2mm*1000um/mm

%                                     rfs{counter}{r}(c) = sqrt(datarun.vision.sta_fits{cell_index(c)}.sd(1)*datarun.vision.sta_fits{cell_index(c)}.sd(2))*2*stixel_size*5*6.5/scale;
                                    
                                end
                                
                            end
                        end
                        
                        
                    end
                catch
                    a=1%                     disp('AHHHH why me no work?')
                end
                
            end
            
            
        end
    end
    
end

%% Process RFs
for i = 1:4
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

figure; plot(rf_cell_type_include(rf_cell_type_include(:,1)>20,1), rf_cell_type_include(rf_cell_type_include(:,1)>20,2), 'ok', 'MarkerFaceColor', 'k')
hold on
xlim = get(gca, 'xlim');
ylim = get(gca,'ylim');
plot([0 max(xlim(2), ylim(2))], [0 max(xlim(2), ylim(2))], '--k')
xlabel('ON parasol RF diameter (\mum)')
ylabel('OFF parasol RF diameter (\mum)')
title(['ON parasols are larger than OFF parasols on average by ', num2str(round(pec_inc_parasol*1000)/10-100), '%'])


figure; plot(rf_cell_type_include(rf_cell_type_include(:,1)>20,3), rf_cell_type_include(rf_cell_type_include(:,1)>20,4), 'ok', 'MarkerFaceColor', 'k')
hold on
xlim = get(gca, 'xlim');
ylim = get(gca,'ylim');
plot([0 max(xlim(2), ylim(2))], [0 max(xlim(2), ylim(2))], '--k')
xlabel('ON midget RF diameter (\mum)')
ylabel('OFF midget RF diameter (\mum)')
title(['ON midgets are larger than OFF midgets on average by ', num2str(round(pec_inc_midget*1000)/10-100), '%'])


for i = 1:size(good_data_sets(:,6),1)
    if isnumeric(good_data_sets{i,6}) == 0
        ecc(i) = 0;
    else
        ecc(i) = good_data_sets{i,6};
    end
    
end

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
plot(ecc(rf_cell_type(:,1)>20), rf_cell_type(rf_cell_type(:,1)>20,1), 'ko', 'MarkerFaceColor', 'w')
plot(ecc(rf_cell_type(:,2)>20), rf_cell_type(rf_cell_type(:,2)>20,2), 'ko', 'MarkerFaceColor', 'k')

xlabel('Temporal equivalent eccentricity (mm)')
ylabel('RF diameter (\mum)')
title('Parasols')
legend('ON parasol', 'OFF parasol');

for i =1:size(good_data_sets,1)
    if rf_cell_type(i,3)<80 && ecc(i) < 8
        rf_cell_type(i,3:4) = nan;
    end
end

figure;
% plot(ecc_include(~isnan(ecc_include)), rf_cell_type_include(~isnan(ecc_include),1), 'ko', 'MarkerFaceColor', 'k')
hold on
plot(ecc(rf_cell_type(:,3)>20), rf_cell_type(rf_cell_type(:,3)>20,3), 'ko', 'MarkerFaceColor', 'w')
plot(ecc(rf_cell_type(:,4)>20), rf_cell_type(rf_cell_type(:,4)>20,4), 'ko', 'MarkerFaceColor', 'k')

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


