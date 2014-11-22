f1=@(x) 1.7 .*x.*(x>=0);
f2=@(x) 1 .*x.*(x>=0);
N= @(x) 1./(1+exp(-5*x)) -0.5;


xGrid = [-0.5:0.01:0.5];
yGrid = [-0.5:0.01:0.5];

xWts=(1/sqrt(2*pi))*exp(-0.5*(xGrid).^2/(0.15)^2);
yWts=xWts;

X_pts = zeros(length(xGrid),length(yGrid));
Y_pts = zeros(length(xGrid),length(yGrid));
Z_pts = zeros(length(xGrid),length(yGrid));
Wt_pts=zeros(length(xGrid),length(yGrid));
inputs=[];
for iptX = 1:length(xGrid)
    for iptY =1:length(yGrid)
        x=xGrid(iptX);
        y=yGrid(iptY);
        
        X_pts(iptX,iptY)=x;
        Y_pts(iptX,iptY)=y;
        Z_pts(iptX,iptY)=N(f1(x)+f2(y));
        Wt_pts(iptX,iptY)=xWts(iptX)*yWts(iptY);
    
    end
end

figure;
subplot(2,2,1);
plot(xGrid,f1(xGrid),'r');
hold on
plotyy(yGrid,f2(yGrid),xGrid,xWts);

subplot(2,2,3);
plot(xGrid,N(xGrid));

subplot(2,2,4);
contourf(X_pts,Y_pts,Z_pts);
colorbar
hold on;
STAx=mean(X_pts(:).*Z_pts(:).*Wt_pts(:));
STAy=mean(Y_pts(:).*Z_pts(:).*Wt_pts(:));
scale=0.5/max(abs(STAx),abs(STAy));
plot([0,scale*STAx,-scale*STAx],[0,scale*STAy,-scale*STAy],'b')
hold on
plot([0,scale*STAx,-scale*STAx],[0,scale*STAy,-scale*STAy],'*')

hold on
len=length(xGrid);

sampleBool=logical(zeros(len,len));
sampleBool(round(len/2),round(len))=1;
sampleBool(round(len/2),round(1.75*len/2))=1;
sampleBool(round(len/2),round(1.5*len/2))=1;
sampleBool(round(len/2),round(1.25*len/2))=1;

sampleBool(round(len),round(len/2))=1;
sampleBool(round(1.75*len/2),round(len/2))=1;
sampleBool(round(1.5*len/2),round(len/2))=1;
sampleBool(round(1.25*len/2),round(len/2))=1;



STAx=mean(mean(X_pts(sampleBool).*Z_pts(sampleBool)))
STAy=mean(mean(Y_pts(sampleBool).*Z_pts(sampleBool)))
scale=0.5/max(abs(STAx),abs(STAy));
plot([0,scale*STAx,-scale*STAx],[0,scale*STAy,-scale*STAy],'r')
hold on
plot([0,scale*STAx,-scale*STAx],[0,scale*STAy,-scale*STAy],'*')



hold on;
thresh=1.9
sampleBool=logical(zeros(len,len));
sampleBool(round(thresh*len/2):end,:)=1;
STAx=mean(mean(X_pts(sampleBool).*Z_pts(sampleBool).*Wt_pts(sampleBool)))

sampleBool=logical(zeros(len,len));
sampleBool(:,round(thresh*len/2):end)=1;
STAy=mean(mean(Y_pts(sampleBool).*Z_pts(sampleBool).*Wt_pts(sampleBool)))
scale=0.5/max(abs(STAx),abs(STAy));
plot([0,scale*STAx,-scale*STAx],[0,scale*STAy,-scale*STAy],'m')
hold on
plot([0,scale*STAx,-scale*STAx],[0,scale*STAy,-scale*STAy],'*')
axis image
% 
% subplot(2,2,2);

subplot(2,2,2);
imagesc(xGrid,yGrid,Wt_pts)
axis image
colorbar


figure;
imagesc(xGrid,yGrid,sampleBool)

