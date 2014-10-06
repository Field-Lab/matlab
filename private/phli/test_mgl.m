%% Pick B&W or RGB
% dims = [600 600];
% dims = [600 600 3];


%% Pick standard, inverted loops, or max optimized
privateFunc = @mglPrivateCreateTexture;
% privateFunc = @mglPrivateCreateTextureGLFormatOpt;

% Needed to get params loaded
if isempty(mglGetParam('frameRate')), mglOpen; mglClose; end

fieldwidth  = 600;
fieldheight = 600;
stixelwidth  = 1;
stixelheight = 1;

x = fieldwidth/stixelwidth;
y = fieldheight/stixelheight;
texx = fieldwidth*mglGetParam('xPixelsToDevice');
texy = fieldheight*mglGetParam('yPixelsToDevice');

t = 60*2;
interval = 1;
framechanges = mglGetParam('frameRate')*t/interval;
framesperttl = 100;
nttls = floor(framechanges/framesperttl);
remainder = mod(framechanges,framesperttl);

predicted_ttl_interval = 1 / mglGetParam('frameRate') * interval * framesperttl;


%%
% clear texture A
% A(:,:,1) = [255 0; 0 0];
% A(:,:,2) = [0 255; 0 0];
% A(:,:,3) = [0 0; 255 0];
% A = matrix_scaled_up(A, [], struct('scale_x', 600, 'scale_y', 300));

% for i = 1:60
%     A{i} = randi([0 255], [600 600 3]);
% end

% A = zeros(600, 600, 4);


%% Matlab frame generation
ttls = zeros(1,nttls+1);
privateFunc = @mglPrivateCreateTexture;
mglOpen; 
try
    tic;
    for t = 1:nttls
        ttls(t) = toc;
        
        for i = 1:framesperttl
            mglClearScreen;
            A = randi([0 1], [3 x y], 'uint8').*245 + 5;
            A(4,:,:) = ones([x y], 'uint8').*255;
            texture = mglCreateTextureDebug(A, privateFunc, [], 0, {'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});
            for i = 1:interval
                mglBltTexture(texture,[0 0 texx texy]);
                mglFlush;
            end
            mglDeleteTexture(texture);
        end
    end
    ttls(end) = toc;
    for i = 1:remainder
        mglClearScreen;
        A = randi([0 1], [3 x y], 'uint8').*245 + 5;
        A(4,:,:) = ones([x y], 'uint8').*255;
        texture = mglCreateTextureDebug(A, privateFunc, [], 0, {'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});
        for i = 1:interval
            mglBltTexture(texture,[0 0 texx texy]);
            mglFlush;
        end
        mglDeleteTexture(texture);
    end
    
    toc
catch err
    err
    mglClose
end
mglClose

ttl_intervals = diff(ttls);
frames_dropped = (ttl_intervals - predicted_ttl_interval) ./ (interval / mglGetParam('frameRate'));
round(frames_dropped(:))
hist(frames_dropped, round(min(frames_dropped)):round(max(frames_dropped)));


%% Vision frame generation
ttls = zeros(1,nttls+1);
rng_enum = get_enum('edu.ucsc.neurobiology.vision.stimulus.FrameGenerator$RandomNumberGenerator', 'JAVA_RANDOM_V2');
color_type =  get_enum('edu.ucsc.neurobiology.vision.stimulus.FrameGenerator$ColorType', 'INDEPENDENT');
fg = edu.ucsc.neurobiology.vision.stimulus.BinaryFrameGenerator(x,y,0.48,rng_enum,11111,color_type);

mglOpen;
try
    tic;
    for t = 1:nttls
        ttls(t) = toc;
        
        for i = 1:framesperttl
            mglClearScreen;
            A = reshape(uint8(fg.nextMGLFrame()), [4 x y]);
            texture = mglCreateTextureDebug(A, privateFunc, [], 0, {'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});
            for i = 1:interval
                mglBltTexture(texture,[0 0 texx texy]);
                mglFlush;
            end
            mglDeleteTexture(texture);
        end
    end
    ttls(end) = toc;
    for i = 1:remainder
        mglClearScreen;
        A = reshape(uint8(fg.nextMGLFrame()), [4 x y]);
        texture = mglCreateTextureDebug(A, privateFunc, [], 0, {'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});
        for i = 1:interval
            mglBltTexture(texture,[0 0 texx texy]);
            mglFlush;
        end
        mglDeleteTexture(texture);
    end
    
    toc
catch err
    err
    mglClose
end
mglClose

ttl_intervals = diff(ttls);
frames_dropped = (ttl_intervals - predicted_ttl_interval) ./ (interval / mglGetParam('frameRate'));
round(frames_dropped(:))
hist(frames_dropped, round(min(frames_dropped)):round(max(frames_dropped)));


%% Asyncronous Vision frame generation
ttls = zeros(1,nttls+1);

rng_enum = get_enum('edu.ucsc.neurobiology.vision.stimulus.FrameGenerator$RandomNumberGenerator', 'JAVA_RANDOM_V2');
color_type =  get_enum('edu.ucsc.neurobiology.vision.stimulus.FrameGenerator$ColorType', 'INDEPENDENT');
fg = edu.ucsc.neurobiology.vision.stimulus.BinaryFrameGenerator(x,y,0.48,rng_enum,11111,color_type);

bufferLength = 120;
mglfa = edu.ucsc.neurobiology.vision.matlab.Matlab.getMGLFrameAccumulator(fg, bufferLength);
th = java.lang.Thread(mglfa);
th.start();

mglOpen;
try
    tic;
    for t = 1:nttls
        ttls(t) = toc;
        
        for f = 1:framesperttl
            mglClearScreen;
            A = reshape(uint8(mglfa.take()), [4 x y]);
            texture = mglCreateTextureDebug(A, privateFunc, [], 0, {'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});
            for i = 1:interval
                mglBltTexture(texture,[0 0 texx texy]);
                mglFlush;
            end
            mglDeleteTexture(texture);
        end
    end
    ttls(end) = toc;
    mglClearScreen;
    A = reshape(uint8(mglfa.take()), [4 x y]);
    texture = mglCreateTextureDebug(A, privateFunc, [], 0, {'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});
    for i = 1:interval
        mglBltTexture(texture,[0 0 texx texy]);
        mglFlush;
    end
    mglDeleteTexture(texture);
    
    toc
catch err
    err
    mglClose
    th.interrupt();
end
mglClose
th.interrupt();
clear mglfa th

ttl_intervals = diff(ttls);
frames_dropped = (ttl_intervals - predicted_ttl_interval) ./ (interval / mglGetParam('frameRate'));
round(frames_dropped(:))
hist(frames_dropped, round(min(frames_dropped)):round(max(frames_dropped)));


%% Manually check scale up
rng_enum = get_enum('edu.ucsc.neurobiology.vision.stimulus.FrameGenerator$RandomNumberGenerator', 'JAVA_RANDOM_V2');
color_type =  get_enum('edu.ucsc.neurobiology.vision.stimulus.FrameGenerator$ColorType', 'INDEPENDENT');
fg = edu.ucsc.neurobiology.vision.stimulus.BinaryFrameGenerator(x,y,0.48,rng_enum,11111,color_type);
A = reshape(uint8(fg.nextMGLFrame()), [4 x y]);
A = permute(A(1:3,:,:), [3 2 1]);
imshow(A); axis xy