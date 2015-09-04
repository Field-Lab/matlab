
load('/Users/alexth/Desktop/backups/backup_paper_dropbox/20141007/3_changing responses/classification/my_data/response_classification_white.mat')
ch_whiteOFF=ch;
load('/Users/alexth/Desktop/backups/backup_paper_dropbox/20141007/3_changing responses/classification/my_data/response_classification_black.mat')
ch_blackOFF=ch;

load('/Users/alexth/Desktop/backups/backup_paper_dropbox/20141007/3_changing responses/classification/my_data/response_classification_white_ONcells.mat')
ch_whiteON=ch;
load('/Users/alexth/Desktop/backups/backup_paper_dropbox/20141007/3_changing responses/classification/my_data/response_classification_black_ONcells.mat')
ch_blackON=ch;

clear ch

for i=4:7
    same_resp_type(i-3) = numel(find(sum(abs(squeeze(ch_whiteOFF(2:3,i,:))-squeeze(ch_whiteOFF(2:3,i+1,:))))==0));  
    dif_resp_type(i-3) = numel(find(sum(abs(squeeze(ch_whiteOFF(2:3,i,:))-squeeze(ch_whiteOFF(2:3,i+1,:))))==1)); 
end

frac_dif_OFF = dif_resp_type./(dif_resp_type + same_resp_type);


for i=4:7
    same_resp_type(i-3) = numel(find(sum(abs(squeeze(ch_blackON(2:3,i,:))-squeeze(ch_blackON(2:3,i+1,:))))==0));  
    dif_resp_type(i-3) = numel(find(sum(abs(squeeze(ch_blackON(2:3,i,:))-squeeze(ch_blackON(2:3,i+1,:))))==1)); 
end

frac_dif_ON = dif_resp_type./(dif_resp_type + same_resp_type);


for i=4:7
    same_resp_type(i-3) = numel(find(sum(abs(squeeze(ch_whiteOFF(2:3,i,:))-squeeze(ch_whiteOFF(2:3,i+1,:))))==0)) + ...
        numel(find(sum(abs(squeeze(ch_blackON(2:3,i,:))-squeeze(ch_blackON(2:3,i+1,:))))==0));  
    dif_resp_type(i-3) = numel(find(sum(abs(squeeze(ch_whiteOFF(2:3,i,:))-squeeze(ch_whiteOFF(2:3,i+1,:))))==1)) + ...
        numel(find(sum(abs(squeeze(ch_blackON(2:3,i,:))-squeeze(ch_blackON(2:3,i+1,:))))==1)); 
end

frac_dif = dif_resp_type./(dif_resp_type + same_resp_type);

figure
for i=1:4
    subplot(2,2,i)
    pie([frac_dif(i), 1-frac_dif(i)], [1 0])
end




