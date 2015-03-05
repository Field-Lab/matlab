function refit(currentLFs, fitParams,currentUnit)

global fitResults

ginput(1);
trial=fitParams.subplots(fitParams.subplots(:,1)==gca,2);


if fitParams.SpikeCount(trial)>5
    borders=fitParams.leftBorder(trial):fitParams.rightBorder(trial);
    
    currentFit=currentLFs(borders,trial);
    currentFit=currentFit*fitParams.currentPolarity;
    
    fit_res=fit(borders',currentFit,'gauss1',...
        'Lower',fitParams.limits(1,1:3),'Upper',fitParams.limits(2,1:3),...
        'Startpoint',[fitParams.amp(trial), fitParams.ind(trial), 30]);
    resFit(1)=fitParams.currentPolarity*fit_res.a1;
    resFit(2)=fit_res.b1;
    resFit(3)=fit_res.c1;
end

fitResults(:,trial,currentUnit)=resFit;

subplot(fitParams.rows,fitParams.cols,find(fitParams.subplots(:,1)==gca))
hold on
x=(fitParams.leftBorder(trial):fitParams.rightBorder(trial))';
y=resFit(1)*exp(-((x-resFit(2))/resFit(3)).^2);
plot(x,y,'r','LineWidth',2)
a=int2str(round(resFit(:))');
a(regexp(a,'  '))='';
a(regexp(a,'  '))='';
text(200,0.8*fitParams.amp(trial),a,'FontSize',10);

display('DONE')


