mov_orignial = mov_orig2;
mov_modify_new=mov_new2;

mov_len=size(mov_orignial,3);
% Variance for whole movie
sd_orig=[];
sd_null=[];

    mask_full=logical(repmat(mask,[1,1,mov_len]));
    mov_orig_mask=mov_orignial(mask_full);
    mov_null_mask=mov_modify_new(mask_full);
    
    sd_orig=[sd_orig;sqrt(var(mov_orig_mask(:)))];
    sd_null=[sd_null;sqrt(var(mov_null_mask(:)))];
% 
% % 
% % [N1,X1]=hist(sd_orig,20);
% % [N2,X2]=hist(sd_null,20);
% % figure;
% % bar(X1,N1,'r');
% % hold on;
% % bar(X2,N2,'b');
% % legend('Original','Null');
% figure('Color','w');
% plot(sd_orig,sd_null,'*');
% xlabel('Original');
% ylabel('Null');
% hold on;
% plot([0,min(sd_orig),max(sd_orig)],[0,min(sd_orig),max(sd_orig)],'r');
% legend('Orig v/s Null','45d line');
% title('Overall variance');
