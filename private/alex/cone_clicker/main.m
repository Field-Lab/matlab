global myPicture cone_list
cone_list=[];
myPicture=imread('/Volumes/Data/Auxiliary/2014-04-10-2/spot/Image1_edited.jpg');
myPicture=myPicture(:,:,1);

hf=figure;
set(hf,'Toolbar','figure')
colormap gray

sbh=subplot(1,1,1);
imagesc(myPicture)

% crop and find all electrodes positions, scale figure

hcrop=uicontrol('style','pushbutton','Units','normalized','position',[0.78 0.945 0.05 0.03],'string','crop','fontsize',16,...
    'callback','crop(sbh)');

hclick=uicontrol('style','pushbutton','Units','normalized','position',[0.6 0.945 0.05 0.03],'string','cones','fontsize',16,...
    'callback','click_cones(sbh,hclick)');




%% ui controls
% get cell from fit
hPickType=uicontrol('style','edit','Units','normalized','position',[0.78 0.945 0.05 0.03],'string','4','fontsize',16,...
    'callback','plot_fits(myFits,datarunA,{str2num(get(hPickType,''String''))}, 9,''r'',true, corX,corY);');
hPickTypeText=uicontrol('style','text','Units','normalized','position',[0.78 0.975 0.05 0.02],'string','cell type','fontsize',12,'background',get(gcf,'color'));


hEnterCell=uicontrol('style','edit','Units','normalized','position',[0.85 0.945 0.05 0.03],'string','0','fontsize',16,...
    'callback','prepare_sta(datarunA,mapHandle,coneHandle,staHandle,hEnterCell,hEnterCell2);');
hEnterCellText=uicontrol('style','text','Units','normalized','position',[0.85 0.975 0.05 0.02],'string','cell ID 1','fontsize',12,'background',get(gcf,'color'));


hEnterCell2=uicontrol('style','edit','Units','normalized','position',[0.91 0.945 0.05 0.03],'string','0','fontsize',16,...
    'callback','prepare_sta(datarunA,mapHandle,coneHandle,staHandle,hEnterCell,hEnterCell2);');
hEnterCell2Text=uicontrol('style','text','Units','normalized','position',[0.91 0.975 0.05 0.02],'string','cell ID 2','fontsize',12,'background',get(gcf,'color'));



% color sliders
hr=uicontrol('style','slider','Units','normalized','position',[0.05 0.7 0.01 0.2],'min',0,'max',2,'value',1,...
    'SliderStep',[0.1,0.2],'callback','redraw(myMap,mapHandle,get(hr,''Value''),get(hg,''Value''),get(hb,''Value''))');
hrText=uicontrol('style','text','Units','normalized','position',[0.07 0.8 0.02 0.02],'string','R','fontsize',12,'background',get(gcf,'color'));


hg=uicontrol('style','slider','Units','normalized','position',[0.05 0.45 0.01 0.2],'min',0,'max',2,'value',1,...
    'SliderStep',[0.1,0.2],'callback','redraw(myMap,mapHandle,get(hr,''Value''),get(hg,''Value''),get(hb,''Value''))');
hgText=uicontrol('style','text','Units','normalized','position',[0.07 0.55 0.02 0.02],'string','G','fontsize',12,'background',get(gcf,'color'));

hb=uicontrol('style','slider','Units','normalized','position',[0.05 0.2 0.01 0.2],'min',0,'max',2,'value',1,...
    'SliderStep',[0.1,0.2],'callback','redraw(myMap,mapHandle,get(hr,''Value''),get(hg,''Value''),get(hb,''Value''))');
hbText=uicontrol('style','text','Units','normalized','position',[0.07 0.3 0.02 0.02],'string','B','fontsize',12,'background',get(gcf,'color'));



