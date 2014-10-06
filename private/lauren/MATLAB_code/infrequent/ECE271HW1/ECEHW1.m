clear all
close all

colormap(gray)

zigzag = [0   1   5   6  14  15  27  28;
2   4   7  13  16  26  29  42;
3   8  12  17  25  30  41  43;
9  11  18  24  31  40  44  53;
10  19  23  32  39  45  52  54;
20  22  33  38  46  51  55  60;
21  34  37  47  50  56  59  61;
35  36  48  49  57  58  62  63];


pathToCheetah = 'ECE271HW1/';

load([pathToCheetah 'TrainingSamplesDCT_8.mat'])

cheetah = imread([pathToCheetah 'cheetah.bmp'],'bmp');
cheetahMask = imread([pathToCheetah 'cheetah_mask.bmp'],'bmp');

cheetah = im2double(cheetah);
cheetahMask = im2double(cheetahMask);

imagesc(cheetah);

histCenters = 1:63;


%% identification of 2nd maximum dct coefficients for training data

BGTrainingScalars = zeros(1, size(TrainsampleDCT_BG,1));
FGTrainingScalars = zeros(1, size(TrainsampleDCT_FG,1));

for i = 1:size(TrainsampleDCT_BG,1)
    chunkWO1 = TrainsampleDCT_BG(i,:);
    chunkWO1(1) = 0;
    BGTrainingScalars(i) = find(abs(chunkWO1) == max(abs(chunkWO1)))-1;
end

for i = 1:size(TrainsampleDCT_FG,1)
    chunkWO1 = TrainsampleDCT_FG(i,:);
    chunkWO1(1) = 0;
    FGTrainingScalars(i) = find(abs(chunkWO1) == max(abs(chunkWO1)))-1;
end


%% transformation and identifications of 2nd maximum dct coefficients for each chunk of cheetah

imageChunks = cell(size(cheetah,1)-7, size(cheetah,2)-7);
dctChunks = cell(size(cheetah,1)-7, size(cheetah,2)-7);
dctChunkVectors = cell(size(cheetah,1)-7, size(cheetah,2)-7);
dctChunkScalarsArray = zeros(size(cheetah,1)-7, size(cheetah,2)-7);


for i = 1:size(cheetah,1)-7
    for j = 1:size(cheetah,2)-7
        imageChunks{i,j} = cheetah(i:i+7, j:j+7);
        dctChunks{i,j} = dct2(imageChunks{i,j});
        chunkWO1 = dctChunks{i,j};
        chunkWO1(1,1) = 0;
        dctChunkScalarsArray(i,j) = zigzag(abs(chunkWO1) == max(max(abs(chunkWO1))));
%        dctChunkVectors{i,j} = zeros(1,64);
%         for k = 1:8
%             for l = 1:8
%                 dctChunkVectors(zigzag(k,l)+1) = dctChunks(k,l);
%             end
%         end
    end
end

dctChunkScalars = reshape(dctChunkScalarsArray,1,size(dctChunkScalarsArray,1)*size(dctChunkScalarsArray,2));


%% Calculating P(x|cheetah) and P(x|background)

BGConProbs = hist(BGTrainingScalars, histCenters)/length(BGTrainingScalars);
FGConProbs = hist(FGTrainingScalars, histCenters)/length(FGTrainingScalars);
cheetahHistogram = hist(dctChunkScalars, histCenters);

%% Calculating P(cheetah) and P(background)

%cheetahProb = sum(sum(cheetahMask))/(size(cheetahMask,1)*size(cheetahMask,2));
cheetahProb = 250/(1053+250);
bgProb = 1 - cheetahProb;

%% plotting histograms

figure
subplot(2,1,1)
bar(BGConProbs);
xlim([0 64])
ylim([0 0.5])
title('P_{X|Y}(x|grass)')

subplot(2,1,2)
bar(FGConProbs);
xlim([0 64])
ylim([0 0.5])
title('P_{X|Y}(x|cheetah)')

% subplot(3,1,3)
% hist(dctChunkScalars, histCenters);
% xlim([1 64])


%% determining best classification scheme using BDR

i_x = zeros(1,63);

for i = 1:63
    pick0 = log(BGConProbs(i)) + log(bgProb);
    pick1 = log(FGConProbs(i)) + log(cheetahProb);
    if pick0 < pick1
        i_x(i) = 1;
    else
        i_x(i) = 0;
    end
end

%% classifying cheetah image

newMask = zeros(size(dctChunkScalarsArray));
for i = 1:size(dctChunkScalarsArray,1)
    for j = 1:size(dctChunkScalarsArray,2)
        newMask(i,j) = i_x(dctChunkScalarsArray(i,j));
    end
end

figure
colormap(gray(255))
imagesc(newMask)
title('calculated mask (array A)')

%% determining probability of error R*

%calculating P(x) for for each i
% 
% p_x0 = 0;
% p_x1 = 0;
% for i = 1:63
%     if i_x == 1
%         p_x1 = p_x1 + BGConProbs(i)*bgProb;
%     else
%         p_x0 = p_x0 + FGConProbs(i)*cheetahProb;
%     end
% end
% 
% %calculating P(i|x)
% p_cheetahx = (FGConProbs./p_x)*cheetahProb;
% p_grassx = (BGConProbs./p_x)*bgProb;
% 
% %calculating risk
% risk = 0;
% for i = 1:63
%     if i_x(i) == 1; %cheetah would chosen for this value of x
%         risk = risk + p_x(i)*p_grassx(i);
%     else
%         risk = risk + p_x(i)*p_cheetahx(i);
%     end
% end

%% calculating error
error = zeros(7,7);

for i = 1:7
    for j = 1:7
        errorMask = abs(cheetahMask(i:i+247, j:j+262) - newMask);
        error(i,j) = sum(sum(errorMask));
    end
end

figure
surf(1:7,1:7,error) %shows dependence of error on shifting the frame of the calculated mask relative to the true mask

%alignment that gives lowest value of error
errorMask = abs(cheetahMask(7:7+247, 1:1+262) - newMask);
figure
imagesc(errorMask)