%


% run once to set variable 'match_list'


if exist('match_list','var')
    
    % identify all potential matches
    match_matrix = [];
    for mm = 1:length(match_list)
        
        axon_id = match_list{mm}{1};
        
        % loop through the matching EIs (will not loop if empty)
        for nn = 1:length(match_list{mm}{2})
            % add each EI that's a good match
            match_matrix = [match_matrix; axon_id match_list{mm}{2}{nn} 1];
        end
        for nn = 1:length(match_list{mm}{3})
            % add each EI that's an ok match
            match_matrix = [match_matrix; axon_id match_list{mm}{3}{nn} 0];
        end
        
    end
    
    
    
    
    
    switch 1
        case 1  %  SLIDER TO SCROLL THROUGH CANDIDATE MATCHES

            % set up figure
            figure(1);clf;
            start_index=1; index_min=1; index_max=size(match_matrix,1);
            slider = make_loop_slider_list(start_index,index_min,index_max);

            while 1
                % create plot
                k = round(get(slider,'Value'));
                subplot('Position',[.05 .1 .90 .85]);cla;

                % get ei, axon ids
                ei_id = match_matrix(k,2);
                axon_ids = match_matrix(k,1);

                % put quality of match in title
                switch match_matrix(k,3)
                    case 1; title_text = 'good match, ';
                    case 0; title_text = 'ok match, ';
                end

                % show what's coming
                title(sprintf('...%saxon %d, EI %d...',title_text,axon_ids,ei_id));drawnow

                % plot EI and axon
                qd_plot_aligned_rgc(datarun,axons,ei_id,axon_ids,'foa',gca,'title',title_text)

                uiwait;
            end


        case 2  % PLOT ALL IN A MATRIX OF SUBPLOTS
            
            % choose match list
            match_list_ = [2086 2181 2191 2197 2296 2311 2312 2314 2418 2431];
            
            % set up axes
            figure(2);clf
            plot_axes = subplot_axes(2,[0 0 1 1],0.1,0.1,3,4);

            for pp = 1:length(match_list_)
                
                % choose ei, axon ID
                ei_id = match_list_(pp);
                axon_ids = 8;
                
                % plot EI and axon
                qd_plot_aligned_rgc(datarun,axons,ei_id,axon_ids,'foa',plot_axes{pp},'title',title_text,...
                    'xlim',[-1000 -350],'ylim',[-500 -120])
                
                %set(plot_axes{pp},'visible','off')
                set(plot_axes{pp},'XTick',[],'ytick',[])
            end
            
            % plot STA fits
            if 0
                datarun = load_sta(datarun,'sync_cell_ids',0,'load_sta',match_list_);
                datarun = load_params(datarun,'load_sta_fits',match_list_,'sync_cell_ids',0,'load_cell_types',0);
                datarun = get_sta_fits_from_vision(datarun,match_list_);
                figure;plot_rf_summaries(datarun,match_list_,'label',1);
            end
    end


else

    switch 1
        case 1

            % 2007-09-18-4/data002-nwpca/data002-nwpca-duplicates/data002
            % excel axons

            %{ axon #, {good matches}, {ok matches}}
            match_list = {...
                { 6, {3601, 3603, 3604}, {3618, 3512} },...
                { 7, {}, {} },...
                { 8, {2296}, {2181} },...
                { 9, {2899}, {2902, 2898, 2911} },...
                { 10, {3020, 3023}, {3017} },...
                { 11, {3258}, {3277} },...
                { 12, {3635, 3619}, {3616, 3618, 3632} },...
                { 14, {2236, 2240}, {2241} },...
                { 15, {2342}, {2346} },...
                { 16, {3203, 3198}, {3196, 3199} },...
                { 19, {3545, 3559}, {3544} },...
                { 21, {}, {3796, 3801} },...
                { 22, {3708}, {3693, 3692, 3698} },...
                { 23, {1801, 1807}, {1743, 1805} },...
                { 23, {1801}, {1804, 1805} },...
                { 24, {1802}, {2992, 2993} },...
                { 26, {}, {3708} },...
                { 27, {}, {3708} },...
                { 31, {1173}, {1291, 1298} },...
                { 32, {}, {1231, 1173} },...
                { 33, {}, {1279} },...
                { 34, {1383}, {1386} },...
                { 35, {}, {1324, 1325} },...
                { 36, {}, {} },...
                { 37, {1201}, {} },...
                { 38, {}, {4311, 4312, 4237} },...
                { 41, {}, {4201, 4202, 4203, 4207} },...
                { 42, {4594, 4791}, {} },...
                { 43, {}, {4597} },...
                { 44, {4400, 4459, 4460}, {4457} },...
                { 45, {4322, 4325}, {} },...
                { 46, {4446}, {} },...
                { 47, {}, {4381, 4384} },...
                { 52, {4759, 4882}, {4876} },...
                { 54, {4744, 4747}, {4743, 4861, 4862, 4865} },...
                { 56, {413}, {} },...
                { 58, {5131, 5013, 5015, 5017}, {5136, 5138} },...
                { 59, {5061, 4942}, {5062} },...
                { 60, {5116, 5123, 5179}, {5117, 5119, 5120} },...
                { 62, {}, {5221} },...
                { 69, {}, {7323} },...
                { 70, {1, 6723}, {8} },...
                { 72, {6481, 6488}, {6487} },...
                { 73, {5557, 5614}, {5556, 5551, 5494, 5498} },...
                { 74, {}, {6122, 6125} },...
                { 75, {}, {} },...
                { 76, {5461, 5463, 5581, 5583}, {5466, 5586, 5587} },...
                { 84, {7126, 7127}, {7128, 7141, 7142, 7143} },...
                { 85, {6543, 6529, 6528}, {6548, 6531, 6526} },...
                { 91, {5836, 5838}, {5826, 5840, 5842} },...
                { 90, {6452, 6572}, {6571} },...
                };


        case 2

            % '2007-09-18-6/data002-nwpca-all/data002-nwpca-all'
            match_list = {...
                { 79, {706, 710, 711}, {766} },...
                { 77, {586, 587}, {590, 592, 593} },...
                { 76, {578}, {572, 573, 577} },...
                { 75, {575}, {} },...
                { 110, {1188}, {1186, 1187, 1189, 1192} },...
                { 109, {1188}, {1186} },...
                { 72, {623, 376, 496, 499, 503}, {501} },...
                { 59, {377, 196}, {201, 198, 202} },...
                { 71, {379, 376}, {} },...
                { 73, {393}, {398} },...
                { 74, {442}, {439} },...
                { 105, {856, 859, 862}, {} },...
                { 104, {}, {1087} },...
                { 102, {}, {4912} },...
                { 111, {1430, 1427, 1428, 1431}, {1429} },...
                { 123, {1536, 1488}, {1487, 1607, 1535} },...
                { 122, {1532, 1533}, {} },...
                { 124, {1727, 1730}, {1849, 1852, 1853} },...
                { 121, {1591, 1594, 1643}, {1711, 1638} },...
                { 108, {}, {} },...
                { 116, {}, {1325, 4488} },...
                { 117, {1503}, {} },...
                { 95, {931, 933, 937}, {935, 1007, 1011} },...
                { 106, {1052, 1114, 1115}, {} },...
                { 94, {998}, {} },...
                { 90, {}, {} },...
                { 93, {}, {} },...
                { 108, {}, {} },...
                { 91, {621}, {} },...
                { 114, {4537, 4598}, {} },...
                { 128, {}, {4414, 4291, 4294, 4298} },...
                { 115, {4593}, {4594} },...
                { 89, {5387, 5789, 5790}, {} },...
                { 69, {123}, {} },...
                { 68, {5552, 5555}, {5551, 5554} },...
                { 101, {5176}, {5061} },...
                { 100, {4997}, {} },...
                { 33, {7487, 7490, 7493, 7367, 7373}, {7473, 7491} },...
                { 44, {7565}, {7593} },...
                { 42, {7233, 7236}, {7096, 7099, 7102, 7111} },...
                { 32, {7141, 7144}, {7143, 7156} },...
                { 29, {7156}, {7171, 7177} },...
                { 31, {7397, 7413, 7416}, {7417} },...
                { 30, {7277}, {7276} },...
                { 57, {6622, 6637}, {} },...
                { 56, {}, {6485, 6497} },...
                { 39, {6320}, {6440, 6317} },...
                { 27, {6707, 6708, 6709}, {} },...
                { 25, {}, {6230, 6232} },...
                { 28, {7037, 7042, 7051}, {7174, 7188, 7021, 7025} },...
                { 40, {}, {6799, 6800} },...
                { 0, {}, {} },...
                { 0, {}, {} },...
                { 0, {}, {} },...
                { 0, {}, {} },...
                { 0, {}, {} },...
                { 0, {}, {} },...
                { 0, {}, {} },...
                { 0, {}, {} },...
                };  % no match found: 19 18 26 41 43 96 97 99 85 86 125



        case 3
            % '2007-01-23-4/data001/data001'
            match_list = {...
                { 1, {}, {4697,4538,4714} },...
                { 2, {2,6841}, {6608} },...
                { 3, {}, {} },...
                { 4, {}, {5587,5528} },...
                };


        case 4
            % '2007-01-23-4/data002/data002'
            match_list = {...
                { 1, {}, {4713,4729,4666} },...
                { 2, {}, {121} },...
                { 3, {}, {} },...
                { 4, {}, {5583} },...
                };



    end


end

