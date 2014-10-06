
% load datarun
if (~exist('a1','var') || ~exist('datarun','var')) && 0
    switch 2

        case 1 % topsirloin
            array_type = 512;
            datarun = load_data('2007-09-18-2/data002-nwpca/data002-nwpca-all/data002');
            datarun = load_java_movie(datarun,'/snle/acquisition/movie-xml/RGB-8-8-0.48-11111.xml');
            a1 = imread('/snle/lab/Experiments/Array/Analysis/2007-09-18-2/images/alive-montage-4x-contrast-flat.tiff');
            a2 = imread('/jacob/snle/data/2007-09-18-2/slidebook/Image 12 - FITC.tif',1);
            f1 = imread('/snle/data/2007-09-18-2/confocal/2009-04-10/2007-09-18-2_Projection043_Merged_ch00.tif');
            % from Series022
            load('/marte/snle/lab/Experiments/Array/Analysis/2007-09-18-2/images/TA1.mat')
            load('/marte/snle/lab/Experiments/Array/Analysis/2007-09-18-2/images/TA1F1.mat')
            load('/marte/snle/lab/Experiments/Array/Analysis/2007-09-18-2/images/TA2.mat')

            ei_id = 4066;  % midget on 272 matches 4066, 4188
            %ei_id = 6841; % cell on 458
            %ei_id = 831;  % cell on 56?
            %ei_id = 5598; % cell on 378?


        case 2 % t-bone
            array_type = 512;
            a1 = imread('/snle/lab/Experiments/Array/Analysis/2007-09-18-4/images/alive-montage-4x-contrast-flat.tiff');
            f1 = imread('/snle/data/2007-09-18-4/confocal/2009-04-13/2007-09-18-4_Projection059_Merged_ch00.tif'); % from Series033
            f2 = imread('/snle/data/2007-09-18-4/confocal/2009-07-06/2007-09-18-4_63x_stitch_01_max.tif');
            f3 = imread('/snle/data/2007-09-18-4/confocal/2009-07-07/2007-09-18-4_63x_stitch_02_relevant/2007-09-18-4_63x_stitch_02_relevant_max.tif');
            f4 = imread('/snle/data/2007-09-18-4/confocal/2009-07-09/2007-09-18-4_63x_stitch_03_max.tif');
            f5 = imread('/snle/data/2007-09-18-4/confocal/2009-07-10/2007-09-18-4_63x_stitch_04_max.tif');
            f6 = imread('/snle/data/2007-09-18-4/confocal/2009-07-10/2007-09-18-4_63x_stitch_05_max.tif');
            
            load('/marte/snle/lab/Experiments/Array/Analysis/2007-09-18-4/images/TA1.mat')
            load('/marte/snle/lab/Experiments/Array/Analysis/2007-09-18-4/images/TA1F1.mat')
            load('/marte/snle/lab/Experiments/Array/Analysis/2007-09-18-4/images/TA1F2.mat')

            axon_spreadsheet_path = '/snle/lab/Experiments/Array/Analysis/2007-09-18-4/images/2007-09-18-4-axons.xls';
            axon_source = 'F1';

            switch 2
                case 1  % data000, scotopic
                    datarun = load_data('2007-09-18-4/data000-nwpca/data000-nwpca-duplicates/data000');
                    ei_id = 2885;  % large cell on 193?
                    ei_id = 3018;  % midget by 202?

                case 2  % data002, photopic
                    datarun = load_data('2007-09-18-4/data002-nwpca/data002-nwpca-duplicates/data002');
                    datarun = load_java_movie(datarun,'/snle/acquisition/movie-xml/RGB-8-8-0.48-11111.xml');
                    ei_id = 1291;axon_ids=31;  % midget on 87?
                    ei_id = 1232;axon_ids=31;  % midget on 87?
                    ei_id = 2886;axon_ids=[];  % midget off array by 177?
                    ei_id = 2226;axon_ids=[];  % cells near 149
                    ei_id = 3545;axon_ids=[];  % midget near 237
                    ei_id = 1279;axon_ids=33;  % ? midget by 86
                    ei_id=1801;axon_ids=23;  % midget by 121
            end


        case 3 % porterhouse
            array_type = 512;
            a1 = imread('/snle/lab/Experiments/Array/Analysis/2007-09-18-6/images/alive-montage-4x-contrast-flat.tiff');
            f1 = imread('/snle/data/2007-09-18-6/confocal/2009-04-08/2007-09-18-6_Projection053_Merged_ch00.tif');
            % from Series028
            load('/marte/snle/lab/Experiments/Array/Analysis/2007-09-18-6/images/TA1.mat')
            load('/marte/snle/lab/Experiments/Array/Analysis/2007-09-18-6/images/TA1F1.mat')

            axon_spreadsheet_path = '/snle/lab/Experiments/Array/Analysis/2007-09-18-6/images/2007-09-18-6-axons-a.xls';
            axon_source = 'F1';

            datarun = load_data('2007-09-18-6/data002-nwpca-all/data002-nwpca-all');
            datarun.stimulus.movie_xml = '/snle/acquisition/movie-xml/RGB-4-48-0.48-11111.xml';


        case 4 % tenderloin
            array_type = 512;
            a1 = imread('/snle/lab/Experiments/Array/Analysis/2007-09-18-7/images/alive-montage-4x-contrast-flat.tiff');
            f1 = imread('');
            % no series which has everything in it, will have to wait for high res scans

            datarun = load_data('2007-09-18-7/data002/data002/data002');
            datarun.stimulus.movie_xml = '/snle/acquisition/movie-xml/BW-4-48-0.48-11111.xml';



        case 5 % groundbeef
            array_type = 512;
            a1 = imread('/snle/lab/Experiments/Array/Analysis/2007-01-23-4/images/alive-10x-montage-flat.tif');
            load('/marte/snle/lab/Experiments/Array/Analysis/2007-01-23-4/images/TA1.mat')

            axon_spreadsheet_path = '/snle/lab/Experiments/Array/Analysis/2007-01-23-4/images/2007-01-23-4-axons.xls';
            axon_source = 'A1';

            switch 2
                case 1 % RGB 16-32
                    datarun = load_data('2007-01-23-4/data001/data001');
                case 2 % drifting gratings
                    datarun = load_data('2007-01-23-4/data002/data002');
            end

    end
    datarun = load_neurons(datarun);

    % get electrode coordinates and array bounds
    switch array_type
        case 512
            % electrode positions (transformed to match array_map)
            ep = electrode_positions(512);

            % region around the array
            array_x = 1100*[-1 1];
            array_y = 500*[-1 1];

            datarun.ei.array_bounds_x = array_x;
            datarun.ei.array_bounds_y = array_y;
    end

    % load axon paths
    if ~isempty(axon_spreadsheet_path)

        % read in spreadsheet
        [junk, junk, axons_xls] = xlsread(axon_spreadsheet_path);

        % load coordinates of the path for each axon
        num_axons = size(axons_xls,2)/2;
        axons = cell(num_axons,1);
        for aa=1:num_axons
            start_col = aa*2-1;
            % determine if contains numeric data
            if isnumeric(axons_xls{2,start_col})
                % load entries
                axons{aa}(:,1) = cell2mat({axons_xls{2:end,start_col}});
                axons{aa}(:,2) = cell2mat({axons_xls{2:end,start_col+1}});
                % remove the trailing NaN entries, which correspond to empty fields
                axons{aa} = axons{aa}(~isnan(axons{aa}(:,1)),:);
            end
        end

        % transform to array space
        axons_orig = axons;
        % make transformation, based on which image was used to draw the axons
        switch axon_source
            case 'F1'
                TA1F1_reverse = cp2tform(base_points,input_points,'lwm'); % these should be from the TA1F1 transformation
                T=maketform('composite',TA1F1_reverse,fliptform(TA1));
            case 'A1'
                T=fliptform(TA1);
        end
        % apply to each axon
        for aa=1:num_axons
            if ~isempty(axons{aa})
                axons{aa} = tforminv(T,axons_orig{aa});
            end
        end
    end

end




% by default, don't show either
ei_id = [];axon_ids = [];


% specify axon and ei
%axon_ids = 33; ei_id = 1276;
ei_scale = 1;


array_x = datarun.ei.array_bounds_x;
array_y = datarun.ei.array_bounds_y;


% initialize title text
title_text = '';

% choose image and transform it
switch 1
    case 1  % A1
        [tform_img,xdata,ydata] = imtransform(a1,TA1,'xdata',array_x,'ydata',array_y);
        title_text = [title_text 'image A1'];
    case 2  % A2
        [tform_img,xdata,ydata] = imtransform(a2,TA2,'xdata',array_x,'ydata',array_y,'XYScale',1/2);
        %[tform_img,xdata,ydata] = imtransform(a2,TA2);
        %tform_img = double(tform_img);
        %tform_img = (tform_img < 1000).*tform_img;
        title_text = [title_text 'image A2'];
    case 3  % F1, with ROI
        roi_x = [1400 2700]; roi_y = [1800 3700];
        [tform_img,xdata,ydata] = imtransform(f1(roi_y(1):roi_y(2),roi_x(1):roi_x(2)),...
            maketform('composite',TA1,TA1F1),'vdata',roi_y,'udata',roi_x);
        title_text = [title_text 'image F1 with ROI'];
    case 4  % F1, no ROI
        [tform_img,xdata,ydata] = imtransform(f1,maketform('composite',TA1,TA1F1),'xdata',array_x,'ydata',array_y);
        title_text = [title_text 'image F1'];
end

% initialize figure
figure(13);clf;subplot('Position',[.05 .02 .95 .98])

% plot image
if 1
    imagesc(tform_img,'xdata',xdata,'ydata',ydata);
else
    %imagesc(permute([.7 .7 .7],[1 3 2]),'xdata',array_x,'ydata',array_y)
end
axis image;axis xy;colormap gray;hold on;
set(gca,'ydir','reverse','xdir','reverse','xlim',array_x,'ylim',array_y)

% plot electrodes
if 0
    plot(ep(:,1),ep(:,2),'.','Color',[1 1 0])
    for ee=1:size(ep,1)
        text(ep(ee,1),ep(ee,2),num2str(ee),'Color',[.5 .5 0],'FontSize',10,...
            'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    end
end

% plot axon path
for axon_id = axon_ids
    title_text = [title_text ', axon'];
    if ~isempty(axons{axon_id})
        plot(axons{axon_id}(:,1),axons{axon_id}(:,2),'Color','r')
        title_text = [title_text sprintf(' %d',axon_id)];
    end
end

% plot EI
if ~isempty(ei_id) && 1
    datarun = load_ei(datarun,ei_id);
    plot_ei(datarun,ei_id,0,'alpha',0,'pretty_axes',0,'pos_color',[.2 .2 1],'max_scale',45,'scale',ei_scale,'cutoff',-1)
    %plot_ei_(new_ei{1},datarun.ei.position,0,'alpha',0,'pretty_axes',0,'pos_color',[1 0 0])
    title_text = [title_text sprintf(', EI %d',ei_id)];
end


title(title_text)


% plot all axon paths
if 1
    figure(101);clf;imagesc(tform_img,'xdata',xdata,'ydata',ydata);;axis image; colormap gray; hold on;
    axis image;axis xy;colormap gray;hold on;
    set(gca,'ydir','reverse','xdir','reverse','xlim',array_x,'ylim',array_y)
    for aa=1:length(axons)
        if ~isempty(axons{aa})
            col = (rand(3,1)+0.3)/1.3;
            plot(axons{aa}(:,1),axons{aa}(:,2),'Color',col,'LineWidth',0.2)
            text(axons{aa}(1,1),axons{aa}(1,2),sprintf('axon %d',aa),'Color',col,'FontSize',5,...
                'HorizontalAlignment','Center','VerticalAlignment','Bottom')
        end
    end
    % plot electrodes
    if 1
        ep = datarun.ei.position;
        plot(ep(:,1),ep(:,2),'.','Color',[1 1 0])
        for ee=1:size(ep,1)
            text(ep(ee,1),ep(ee,2),num2str(ee),'Color',[.5 .5 0],'FontSize',3,...
                'HorizontalAlignment','Center','VerticalAlignment','Bottom')
        end
    end
end

% compute and plot STA
if 0
    % compute movie
    if ~exist('vision_movie','var')
        vision_movie = compute_vision_movie(datarun.stimulus.movie_xml,datarun.triggers);end
    % compute STA
    sta = compute_sta_(datarun.spikes{get_cell_indices(datarun,cell_id)}(proj(:,1)>20 & proj(:,2)<0),...
        datarun.triggers, vision_movie);  %(proj(:,1)>20 & proj(:,2)<0)
    % plot
    figure;imagesc(norm_image(rf_from_sta(sta_from_java_sta(sta))));axis image;title(sprintf('RF of cell id %d',cell_id))
end

drawnow


if 0


    % 2007-09-18-4/data002-nwpca/data002-nwpca-duplicates/data002

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
        { 0, {}, {} },...
        };



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




    % '2007-01-23-4/data001/data001'
    match_list = {...
        { 1, {}, {4697,4538,4714} },...
        { 2, {2,6841}, {6608} },...
        { 3, {}, {} },...
        { 4, {}, {5587,5528} },...
        };

    % '2007-01-23-4/data002/data002'
    match_list = {...
        { 1, {}, {4713,4729,4666} },...
        { 2, {}, {121} },...
        { 3, {}, {} },...
        { 4, {}, {5583} },...
        };


end



% old code for aligning image to array image
% % array image
% array_map = imread('/snle/home/gauthier2/Desktop/imaging/methods/alignment graphics/512.tiff');
% % transformation from array image to electrode positions
% switch version('-release')
%     case '2007a'
%         TE=cp2tform([150.5 207.5; 894.5 207.5; 138.5 567.5; 882.5 567.5],ep([249 392 129 512],:),'affine');
%     otherwise
%         TE=cp2tform([150.5 207.5; 894.5 207.5; 138.5 567.5; 882.5 567.5],ep([249 392 129 512],:),'similarity');
% end
