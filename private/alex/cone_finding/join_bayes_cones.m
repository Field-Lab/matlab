
piece = '2012-09-13-2';
run = 'data009';

% define data path
datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',false);
datarun=load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = set_polarities(datarun);
datarun = load_sta(datarun,'load_sta',[]);


datarun_ONp = datarun;
datarun_ONp = load_cones(datarun_ONp,'ONparasol_');

datarun_ONm = datarun;
datarun_ONm = load_cones(datarun_ONm,'ONmidgets_');

datarun_OFFp = datarun;
datarun_OFFp = load_cones(datarun_OFFp,'OFFparasol');

datarun_OFFm = datarun;
datarun_OFFm = load_cones(datarun_OFFm,'OFFmidgets');

datarun_SBC = datarun;
datarun_SBC = load_cones(datarun_SBC,'SBC10_');


cone_data = find_cone_data(datarun_ONp);
myData=detect(cone_data, @(cd) (regexp(cd, 'ONparasol')));
file_prefix=[single_cone_path myData,'/'];
load([file_prefix 'parameters.mat']);
load([file_prefix 'mydll.mat']);
ONparasoldll=mydll;

myData=detect(cone_data, @(cd) (regexp(cd, 'OFFparasol')));
file_prefix=[single_cone_path myData,'/'];
load([file_prefix 'mydll.mat']);
OFFparasoldll=mydll;

myData=detect(cone_data, @(cd) (regexp(cd, 'ONmidgets')));
file_prefix=[single_cone_path myData,'/'];
load([file_prefix 'mydll.mat']);
ONmidgetsdll=mydll;

myData=detect(cone_data, @(cd) (regexp(cd, 'OFFmidgets')));
file_prefix=[single_cone_path myData,'/'];
load([file_prefix 'mydll.mat']);
OFFmidgetsdll=mydll;

myData=detect(cone_data, @(cd) (regexp(cd, 'SBC10')));
file_prefix=[single_cone_path myData,'/'];
load([file_prefix 'mydll.mat']);
SBCdll=mydll;

a=datarun_ONp.cones.centers;
b=datarun_OFFp.cones.centers;
c=datarun_ONm.cones.centers;
d=datarun_OFFm.cones.centers;
e=datarun_SBC.cones.centers;



showMe={'ON parasol', 'OFF parasol','SBC'};
cone_labels=[false false true];

myColorData=zeros(size(OFFparasoldll));
myCones=cell(1,3);
for i=1:3 
    switch showMe{i}
        case 'ON parasol'
            myColorData(:,:,i)=ONparasoldll(:,:,1);
            myCones{i}=a;
        case 'OFF parasol'
            myColorData(:,:,i)=OFFparasoldll(:,:,1);
            myCones{i}=b;
        case 'ON midget'
            myColorData(:,:,i)=ONmidgetsdll(:,:,1);
            myCones{i}=c;
        case 'OFF midget'
            myColorData(:,:,i)=OFFmidgetsdll(:,:,1);
            myCones{i}=d;
        case 'SBC'
            myColorData(:,:,i)=SBCdll(:,:,1);
            myCones{i}=e;
    end
end

figure
imagesc(myColorData,'xdata',[1 datarun.stimulus.field_width]-(1-bcf_params.kernel_spacing),...
    'ydata',[1 datarun.stimulus.field_height]-(1-bcf_params.kernel_spacing))
axis image; hold on
colors=[0.8 0 0; 0 0.8 0; 0 0 0.8];
markers='+xd';

for i=1:3
    if cone_labels(i)
        plot(myCones{i}(:,1),myCones{i}(:,2),markers(i),'color',colors(i,:),'MarkerSize',3)
    end
end










figure
subplot(2,2,1)
hist(min(pdist2(a,b)),0.5:0.5:10)
axis([0 10 0 inf])

subplot(2,2,2)
hist(min(pdist2(c,d)),0.5:0.5:10)
axis([0 10 0 inf])

subplot(2,2,3)
hist(min(pdist2(a,e)),0.5:0.5:10)
axis([0 10 0 inf])

subplot(2,2,4)
hist(min(pdist2(d,e)),0.5:0.5:10)
axis([0 10 0 inf])

figure
clear m
for i=1:length(e)    
    D=pdist2(e(i,:),e);
    m(i)=min(D(D>0));
end
hist(m,5:1:20)
axis([5 20 0 inf])


m=[min(pdist2(a,b)) min(pdist2(b,c)) min(pdist2(c,d)) min(pdist2(d,a))];
figure
hist(m,0.5:0.5:10)
axis([0 10 0 inf])

threshold=2;

allCones=[a; b; c; d; e];
newCones=[];

clear m
for i=1:length(allCones)    
    D=pdist2(allCones(i,:),allCones);
    m(i)=sum(D<threshold & D>0);
    newCones=[newCones; allCones(i,:), m(i)];
    allCones(D<threshold,:)=0;
end
figure
hist(m,0:1:10)
axis([0 10 0 inf])
figure
plot(a(:,1),a(:,2),'xr')
hold on
plot(b(:,1),b(:,2),'xb')
plot(c(:,1),c(:,2),'xm')
plot(d(:,1),d(:,2),'xc')
plot(e(:,1),e(:,2),'xg')
plot(newCones(newCones(:,3)>0,1),newCones(newCones(:,3)>0,2),'.k','markerSize',10)
legend('ON parasol','OFF parasol','ON midget','OFF midget','SBC','shared')





%% MAKE CONE WEIGHTS MATRIX

% size
rf_size = [datarun.stimulus.field_height datarun.stimulus.field_width];

joinedConeCenters=[datarun_SBC.cones.centers; datarun_ONm.cones.centers; datarun_ONp.cones.centers;...
    datarun_OFFm.cones.centers; datarun_OFFp.cones.centers];
totalNcones=size(joinedConeCenters,1);

% get locations
kernel_spec = struct('center',mat2cell(joinedConeCenters,repmat(1,totalNcones,1)));


cone_data = find_cone_data(datarun_OFFp);
myData=detect(cone_data, @(cd) (regexp(cd, 'ONparasol')));
file_prefix=[single_cone_path myData,'/'];
load([file_prefix 'parameters.mat']);


% cone type
a=[datarun_SBC.cones.types; datarun_OFFp.cones.types];
a(a=='U')='C';

for i=1:length(kernel_spec)
    kernel_spec(i).type = a(1);
    kernel_spec(i).radius = bcf_params.kernel_radii(1);
end


a=datarun_ONp.cones.centers;
b=datarun_ONm.cones.centers;
c=datarun_OFFp.cones.centers;
d=datarun_OFFm.cones.centers;
e=datarun_SBC.cones.centers;

figure
plot(a(:,1),a(:,2),'r+')
hold on
plot(b(:,1),b(:,2),'bx')
plot(c(:,1),c(:,2),'go')
plot(d(:,1),d(:,2),'kd')

plot(e(:,1),e(:,2),'.b','markerSize',10)

legend('ON parasol','ON midget','OFF parasol','OFF midget','SBC')


% make matrix
W = make_cone_weights_matrix(rf_size,kernel_spec,bcf_params.kernel_colors);



% get list of cell indices
cell_spec={1 2 3 4 5};
cell_indices = get_cell_indices(datarun,'all');

% note when it started
fprintf('\nComputing cone weights in %d RFs',length(cell_indices));
start_time_regress = clock; 

% initialize
cone_weights = zeros(totalNcones,length(cell_indices));
robust_std_method=1;
% go through list of cells
for cc = 1:length(cell_indices)
    
    fprintf('.')

    % get summary frame
%     temp_frame = length(datarun.stas.time_courses{1}(:,1)) -1;
%     [rf, junk] = get_sta(datarun,datarun.cell_ids(cell_indices(cc)), 'frames', temp_frame);
     rf = get_rf(datarun,datarun.cell_ids(cell_indices(cc)));
    if isempty(rf)
        continue
    end
    
    rf = repmat(rf,1,1,3);

    % reshape for the regression
    rf = reshape(rf,[],1);

    % put in units of SNR
    rf = rf / robust_std(rf, robust_std_method);


    % regress to get the cone weights

    cone_weights(:,cc) = W\rf; 

end

datarun_joined=datarun;
datarun_joined.cones.centers=joinedConeCenters;
datarun_joined.cones.weights=cone_weights;


%%

datarun = load_sta(datarun,'load_sta','all');


bw=datarun_joined.cones.weights;

% SBC
sbc=datarun.cell_types{5}.cell_ids;
ccb=datarun_joined.cones.centers;

for tmp=1:length(sbc)
    
    unit=find(datarun.cell_ids==sbc(tmp),1);
    
    if nnz(aw(:,unit))
            
        figure
        sta=squeeze(datarun.stas.stas{unit}(:,:,1,5));
        colormap gray
        imagesc(sta)
        hold on
        curcones=(bw(:,unit)/max(abs(bw(:,unit)))+1)/2; % normalized to mean=0.5 and max=1
        curconesAbs=abs(bw(:,unit)/max(bw(:,unit))); % normalized to max 1, mean 0, and then absolute
        for i=1:length(ccb)
            
            if curcones(i)>0.5
                col=[curcones(i) 0 0];
                siz=(curconesAbs(i))*40/(mean(curconesAbs)*30)+1;
            else
                col=[0 curcones(i) 0];
                siz=(curconesAbs(i))*40/(mean(curconesAbs)*30)+1;
            end
            plot(ccb(i,1),ccb(i,2),'+','markersize',siz,'color',col);
        end
        
        title(['SBC, i=',int2str(unit)])
    end
end



% OFF midget
sbc=datarun_joined.cell_types{4}.cell_ids;
ccb=datarun_joined.cones.centers;

for tmp=1:length(sbc)
    
    unit=find(datarun.cell_ids==sbc(tmp),1);
    
    if nnz(bw(:,unit))
            
        figure
        sta=squeeze(datarun.stas.stas{unit}(:,:,1,5));
        colormap gray
        imagesc(sta)
        hold on
        curcones=(bw(:,unit)/max(abs(bw(:,unit)))+1)/2; % normalized to mean=0.5 and max=1
        curconesAbs=abs(bw(:,unit)/max(bw(:,unit))); % normalized to max 1, mean 0, and then absolute
        for i=1:length(ccb)
            
            if curcones(i)>0.5
                col=[curcones(i) 0 0];
                siz=(curconesAbs(i))*40/(mean(curconesAbs)*30)+1;
            else
                col=[0 curcones(i) 0];
                siz=(curconesAbs(i))*40/(mean(curconesAbs)*30)+1;
            end
            plot(ccb(i,1),ccb(i,2),'+','markersize',siz,'color',col);
        end
        
        title(['ON parasol, i=',int2str(unit)])
    end
end

