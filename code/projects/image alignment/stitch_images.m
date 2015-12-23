% Stitch together images


% PARAMETERS
switch 15
    case 1
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/Volumes/gauthier/Confocal/2009-07-06/2007-09-18-4/';
        name = @(series,z,ch)sprintf('2007-09-18-4_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-4_63x_stitch_01_z%03d.tif',z);
        z_range = 0:49;
        chs = {[],0,1};
        positions{1} = {72 69 68 65};
        positions{2} = {71 70 67 66};
    case 2
        size_x = 1024;size_y = size_x;
        im_type = 'uint16';
        im_path = '/Volumes/gauthier/Confocal/2009-07-07/2007-09-18-4/';
        name = @(series,z,ch)sprintf('2007-09-18-4_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-4_63x_stitch_02_z%03d.tif',z);
        z_range = 0:52;
        chs = {[],0,1};
        positions{1} = {36 35 32 31 28};
        positions{2} = {37 34 33 30 29};
    case 3
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-4/confocal/2009-07-09/';
        name = @(series,z,ch)sprintf('2007-09-18-4_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-4_63x_stitch_03_z%03d.tif',z);
        z_range = 0:66;
        chs = {[],0,1};
        positions{1} = {27 26 23 22 19};
        positions{2} = {28 25 24 21 20};
    case 4
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-4/confocal/2009-07-10/';
        name = @(series,z,ch)sprintf('2007-09-18-4_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-4_63x_stitch_04_z%03d.tif',z);
        z_range = 0:60;
        chs = {[],0,1};
        positions{1} = {27 26 23 22 19};
        positions{2} = {28 25 24 21 20};
    case 5
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-4/confocal/2009-07-10/';
        name = @(series,z,ch)sprintf('2007-09-18-4_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-4_63x_stitch_05_z%03d.tif',z);
        z_range = 0:59;
        chs = {[],0,1};
        positions{1} = {49 48 45 44 41};
        positions{2} = {50 47 46 43 42};
        
    case 6
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-0/confocal/2009-09-14/';
        name = @(series,z,ch)sprintf('2007-09-18-0_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-0_63x_stitch_02_z%03d.tif',z);
        z_range = 0:69;
        chs = {[],1,0};
        positions{1} = {27 26 23 22 19};
        positions{2} = {28 25 24 21 20};
        
    case 7
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-0/confocal/2009-09-17/';
        name = @(series,z,ch)sprintf('2007-09-18-0_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-0_63x_stitch_03_z%03d.tif',z);
        z_range = 36:59;
        chs = {[],1,0};
        positions{1} = {38 37 34 33 30};
        positions{2} = {39 36 35 32 31};
        
    case 8
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-0/confocal/2009-09-22/';
        name = @(series,z,ch)sprintf('2007-09-18-0_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-0_63x_stitch_04_z%03d.tif',z);
        z_range = 0:61;
        chs = {[],1,0};
        positions{1} = {39 38 35 34 31};
        positions{2} = {40 37 36 33 32};
        
    case 9
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-0/confocal/2009-09-24/';
        name = @(series,z,ch)sprintf('2007-09-18-0_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-0_63x_stitch_05_z%03d.tif',z);
        z_range = 0:62;
        chs = {[],1,0};
        positions{1} = {44 43 40 39 36};
        positions{2} = {45 42 41 38 37};
        
    case 10
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-0/confocal/2009-09-02/';
        name = @(series,z,ch)sprintf('2007-09-18-0_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-0_63x_stitch_01_z%03d.tif',z);
        z_range = 0:40;
        chs = {[],1,0};
        positions{1} = {46 45 42 41 38};
        positions{2} = {47 44 43 40 39};
        
    case 11
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-6/confocal/2010-02-20/';
        name = @(series,z,ch)sprintf('2007-09-18-6_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-6_63x_stitch_01_z%03d.tif',z);
        z_range = 0:47;
        chs = {[],1,0};
        positions{1} = {46 43 42 39};
        positions{2} = {45 44 41 40};
        
    case 12
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-6/confocal/2010-02-20b/';
        name = @(series,z,ch)sprintf('2007-09-18-6_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-6_63x_stitch_02_z%03d.tif',z);
        z_range = 0:73;
        chs = {[],1,0};
        positions{1} = {25 22 21 18};
        positions{2} = {24 23 20 19};
        
    case 13
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-6/confocal/2010-02-21/';
        name = @(series,z,ch)sprintf('2007-09-18-6_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-6_63x_stitch_03_z%03d.tif',z);
        z_range = 0:67;
        chs = {[],1,0};
        positions{1} = {30 27 26 23};
        positions{2} = {29 28 25 24};
        
    case 14
        size_x = 1024;size_y = size_x;
        im_type = 'uint8';
        im_path = '/snle/data/2007-09-18-6/confocal/2010-02-21/';
        name = @(series,z,ch)sprintf('2007-09-18-6_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-6_63x_stitch_04_z%03d.tif',z);
        z_range = 0:56;
        chs = {[],1,0};
        positions{1} = {55 52 51 48};
        positions{2} = {54 53 50 49};
        
    case 15
        
        
        size_x = 512;size_y = 512;
        im_type = 'uint8';
        im_path = '/snle/data/2010-03-31-0/confocal/2010-04-02/';
        name = @(series,z,ch)sprintf('2010-03-31-0_Series%03d_z%03d_ch%02d.tif',series,z,ch);
        new_name =  @(z)sprintf('2007-09-18-6_63x_stitch_04_z%03d.tif',z);
        z_range = 0:56;
        chs = {[],1,0};
        
        pos={};nx=9;ny=13;pt=[nx 1];
        for kk=1:(nx*ny);
            disp(pt);
            pos{pt(2)}{pt(1)}=kk;
            if mod(pt(1),2)==1;
                if pt(2)<ny;pt(2)=pt(2)+1;
                else
                    pt(1)=pt(1)-1;
                end;
            else
                if pt(2)>1;pt(2)=pt(2)-1;
                else
                    pt(1)=pt(1)-1;
                end;
            end;
        end
        positions = pos;
end


% STITCH!

% number of images across and down
num_y = length(positions);
num_x = length(positions{1});

% reconstruct for each z plane and save out individually
for zz = 1:length(z_range)

    % initialize
    big_image = zeros([num_y*size_y num_x*size_x length(chs)],im_type);

    % cycle through each channel
    for cc = 1:length(chs)
        
        % if the channel exists
        if ~isempty(chs{cc})

            % cycle through each series (i.e. each location)
            for yy = 1:num_y
                for xx = 1:num_x
                    % note range
                    x_range = 1+(xx-1)*size_x : xx*size_x;
                    y_range = 1+(yy-1)*size_y : yy*size_y;
                    
                    % load source image
                    source_im = imread([im_path name(positions{yy}{xx},z_range(zz),chs{cc})]);

                    % place in proper location
                    big_image(y_range,x_range,cc) = source_im;
                    
                end
            end
        end
    end

    
    figure(3);clf;image(double(big_image)/max(max(max(double(big_image)))));axis image;drawnow
    
    switch im_type
        case 'uint16' % if 12 bit, reduce image size
            big_image = uint8(big_image/2^4);
        case 'uint8' % if 8 bit, do nothing
        otherwise
            error('unrecognized image type ''%s''',im_type)
    end
    
    % save out with new name
    imwrite(big_image,[im_path new_name(z_range(zz))],'tif')
    
    
end

