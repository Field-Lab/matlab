for i=1:75
%%

%% load desired frame from the binary file
fid = fopen('crest_s2_t2_b_d5.rawMovie');
fseek(fid,0,'bof');
tline = fgetl(fid);
while ~strcmp(tline,'')
  tline = fgetl(fid);
end
fseek(fid,(i-1)*(300*300*3),'cof');
stim_frame = fread(fid,[300*300*3]);
%mat = reshape(stim_frame,[3 300 300]);
%mat = squeeze(mat(1,:,:))' / 256;
%test1 = repmat(mat,[1 1 3]);
%test1raw = stim_frame(:);
stim_frame = stim_frame / 256;
fclose(fid);
%% load desired frame from the mat file
load crest_s2_t2_b_d5.mat
test2 = mat_out(:,:,i);
test2 = repmat(test2,[1 1 3]);
%% compare the two
if 1
clf
subplot(3,1,1);
imagesc(test1); axis image;
subplot(3,1,2);
imagesc(test2); axis image;
subplot(3,1,3);
imagesc(test1(:,:,1)-test2,[-1 1]); axis image;
end

%% compute input in cone space
load Wc.mat
cone_input = Wc'*test1(:);
ind = 1:size(Wc,1);
ind = reshape(ind, 300, 300, 3);
ind = permute(ind, [3 2 1]);
ind = reshape(ind, [], 1);
Wcp = Wc(ind,:);
cone_input2 = Wcp'*test1raw(:);
%% look at stimulus in cone space
clf
imagesc(test2);
hold on
scatter(dat.locations(:,1),dat.locations(:,2),40,cone_input2-0.5,'filled');
pause(0.01);
colormap(redgrayblue);
axis image off;
end
%%
