function su_activation_plot_gamma(fitGMLM,mov_filtered,gamma,su)

Sigma = mov_filtered*mov_filtered'/size(mov_filtered,2);

isu=su(1);
jsu=su(2);

isig = sqrt(fitGMLM.Linear.filter{isu}'*Sigma*fitGMLM.Linear.filter{isu});
jsig = sqrt(fitGMLM.Linear.filter{jsu}'*Sigma*fitGMLM.Linear.filter{jsu});
ijsig = sqrt(fitGMLM.Linear.filter{isu}'*Sigma*fitGMLM.Linear.filter{jsu});
kSig = [isig^2 , ijsig^2;
        ijsig^2 , jsig^2];

x= [0.01:0.03:2];
    
    
cdf_data=zeros(length(x),length(x));

icnt=0;
for ix=x
    icnt=icnt+1;
    jcnt=0;
    for jx=x
        jcnt=jcnt+1;
        xx= [(ix).^(1/gamma) ; (jx).^(1/gamma)];
%         xxy = W*xx;
        cdf_data(icnt,jcnt)=mvncdf(xx,[0;0],kSig);
%        pause
    end
end

%%
cdf_data_diff1=0*cdf_data;
for idim=2:length(x)
    for jdim=1:length(x)
    cdf_data_diff1(idim,jdim) = cdf_data(idim,jdim)-cdf_data(idim-1,jdim);
    end
end


cdf_data_diff2=0*cdf_data;
for idim=2:length(x)
    for jdim=1:length(x)
    cdf_data_diff2(jdim,idim) = cdf_data_diff1(jdim,idim)-cdf_data_diff1(jdim,idim-1);
    end
end

pdf_data = cdf_data_diff2;

%% 
figure;
imagesc(pdf_data);

pdf_data=pdf_data/sum(pdf_data(:));
figure('Color','w');
subplot(2,2,2);
xgrid = repmat(x,[length(x),1]);
ygrid = repmat(x',[1,length(x)]);
contourf(xgrid,ygrid,pdf_data,20);
title('joint pdf su 1 - 2');

subplot(2,2,1);
plot(x,sum(pdf_data,2));
title('pdf su 1');

subplot(2,2,4);
plot(x,sum(pdf_data,1));
title('pdf su 2');

end