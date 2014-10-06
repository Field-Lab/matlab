fieldwidth  = 600;
fieldheight = 600;
stixelwidth  = 15;
stixelheight = stixelwidth;

x = fieldwidth/stixelwidth;
y = fieldheight/stixelheight;

%% BW white noise

% rng_enum = get_enum('edu.ucsc.neurobiology.vision.stimulus.FrameGenerator$RandomNumberGenerator', 'JAVA_RANDOM_V2');
% color_type =  get_enum('edu.ucsc.neurobiology.vision.stimulus.FrameGenerator$ColorType', 'DEPENDENT');
% fg = edu.ucsc.neurobiology.vision.stimulus.BinaryFrameGenerator(x,y,0.48,rng_enum,11111,color_type);
nframes = 60;
F = struct('cdata', [], 'colormap', []);
for i = 1:nframes
    A = randi([0 1], [y x 1], 'uint8').*245 + 5;
    imshow(A, 'InitialMagnification', 2000);
    set(gca, 'Position', [0 0 1 1])
    drawnow();
    pause();
    F(end+1) = getframe();
end
close
F = F(2:end);

filename = sprintf('BW-%d-%dx%d', stixelwidth, fieldwidth, fieldheight);


%% Moving bar

width = 60;
speed = 6;
nframes = fieldwidth/speed + 1;
grey = uint8(255/2);
blank = grey.*ones(fieldwidth, fieldheight, 'uint8');
F = struct('cdata', [], 'colormap', []);
for i = 1:nframes
    pos = (i-1)*speed + 1;
    A = blank;
    A(:,max(pos-width/2,1):min(pos+width/2,fieldwidth)) = 0;
    imshow(A, 'InitialMagnification', 50);
    set(gca, 'Position', [0 0 1 1])
    drawnow();
    F(end+1) = getframe();
end
close
F = F(2:end);

filename = sprintf('movingbar', stixelwidth, fieldwidth, fieldheight);

%% Save
vw = VideoWriter(filename, 'Uncompressed AVI');
vw.open();
vw.writeVideo(F);
vw.close();