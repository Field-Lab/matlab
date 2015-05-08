dx=0.01;
x=[-5:dx:5];
p1x = @(x) 1/(1+abs(x));
y=x;
pxy=[];
icntx=0;
for ix=x;
    icntx=icntx+1;
    icnty=0;
    for iy=y;
       icnty=icnty+1;
       pxy(icntx,icnty)=p1x(ix)*p1x(iy); 
    end
end

 
figure;
contour(x,y,pxy)
CDF(1)=p1x(x(1))*dx;
for icnt=2:length(x)
CDF(icnt)= CDF(icnt-1) + p1x(x(icnt))*dx;   
end
CDF=CDF/max(CDF);
figure('Color','w');
plot(x,CDF);
title('CDF cauchy');

dxinv=0.01;
icnt=1;
CDFinv=[];
for xinv=[0:dxinv:1]
    aa=x(abs(xinv-CDF)<0.005);
CDFinv(icnt) = aa(1);
icnt=icnt+1; 
end


sigma=1;
mu=0;
p2x = @(x) (1/sqrt(2*pi*sigma^2))* (exp(-((x-mu)^2)/(2*sigma^2) ));
CDF2(1)=p2x(x(1))*dx;
for icnt=2:length(x)
CDF2(icnt)= CDF2(icnt-1) + p2x(x(icnt))*dx;   
end
figure('Color','w');
plot(x,CDF2);
title('CDF gaussian');
% 
% y=x;
% pxy=[];
% icntx=0;
% for ix=x;
%     icntx=icntx+1;
%     icnty=0;
%     for iy=y;
%        icnty=icnty+1;
%        pxy(icntx,icnty)=p2x(ix)*p2x(iy); 
%     end
% end
%  
% figure;
% contour(x,y,pxy)

NL = CDF./CDF2;
figure;
plot(x,NL);
title('Non-linearity')

icnt=1;
F2invF1=[];
for ix=x
F1 = CDF2 (ix==x);
    aa=x(abs(F1-CDF)<0.005);
F2invF1(icnt) = aa(1);
icnt=icnt+1;
end

figure('Color','w');
plot(x,F2invF1);
title('Sub-unit Non-linearity');

% Actual SU non linearity
