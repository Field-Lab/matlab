
function su_activation_plot(fitGMLM,mov_filtered)
Sigma = mov_filtered*mov_filtered'/size(mov_filtered,2);


%% make joint cdf
isu=1;
jsu=2;
fi = fitGMLM.Linear.filter{isu};
fj = fitGMLM.Linear.filter{jsu};
f=[fi,fj];
fmax = max(f')';
%f=f.*(f==repmat(fmax,[1,2]));

isig = sqrt(f(:,1)'*Sigma*f(:,1));
jsig = sqrt(f(:,2)'*Sigma*f(:,2));
ijsig = sqrt(f(:,1)'*Sigma*f(:,2));

kSig = [isig^2 , ijsig^2;
        ijsig^2 , jsig^2];
[E,lam] = eig(kSig);
W = (lam^(-0.5)) * E';

ib = fitGMLM.Linear.bias{isu};
jb = fitGMLM.Linear.bias{jsu};

x=[0.01:.05:10];
cdf_data=zeros(length(x),length(x));

icnt=0;
for ix=x
    icnt=icnt+1;
    jcnt=0;
    for jx=x
        jcnt=jcnt+1;
        xx= [(log(ix)-ib) ; (log(jx)-jb)];
%         xxy = W*xx;
        cdf_data(icnt,jcnt)=mvncdf(xx,[0;0],kSig);
%        pause
    end
end

%% make pdf

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
p1 = sum(pdf_data,2);
plot(x,sum(pdf_data,2));
title('pdf su 1');

subplot(2,2,4);
p2 = sum(pdf_data,1);
plot(x,sum(pdf_data,1));
title('pdf su 2');

subplot(2,2,3);
xgrid = repmat(x,[length(x),1]);
ygrid = repmat(x',[1,length(x)]);
pind = repmat(p1',[length(x),1]) .* repmat(p2',[1,length(x)]);
contourf(xgrid,ygrid,pind',20);
title('If sub-units independent');

end
