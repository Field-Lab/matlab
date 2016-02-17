%% Startup paths
% setup cvx
addpath(genpath('/Users/bhaishahster/Downloads/cvx'))
cvx_setup

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');


%% Preprocess the data lauren gave

elecRespData = load('/Volumes/Lab/Users/grosberg/transfer/responseData_2015_11_09_3_data001_002.mat');

nElecs = 512;
cellIDs=[];
for ielec=1:nElecs
cellIDs = [cellIDs;elecRespData.allCells{ielec}];
end
cellIDs = sort(unique(cellIDs),'ascend');
nCells = length(cellIDs);

resp_mat=zeros(nCells,nElecs);
for ielec=1:nElecs
    nc = length(elecRespData(1).allCells{ielec});
    for icell=1:nc
    cell = elecRespData.allCells{ielec}(icell);
    prob = elecRespData.allprobs{ielec}(icell);
    resp_mat(cellIDs==cell,ielec) = prob;
    end
end



WN_datafile = '/Volumes/Analysis/2015-11-09-3/data000/data000';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

%% Select cells and make them fire!

%selected_cells = [4519,4653,4666,4817,4956];
%selected_cells = [3228,3287,3454,3541,3306,3319,2411,3439,3471,3513,3556];
selected_cells = [3454,3306,3319,2411,3439,3471,3513,3556];


% clean cells and make response vector
final_selected_cells=[];resp_vec=zeros(nCells,1);
for icell=1:length(selected_cells)
    
    if(sum(selected_cells(icell)==cellIDs)==1)
        final_selected_cells = [final_selected_cells;selected_cells(icell)];
    end
    resp_vec = resp_vec + double(selected_cells(icell)==cellIDs);
end


cvx_begin
    variable current(512)
    subject to
        RC = resp_mat*current
        current>=0
        current<=1
        sum(current)<=3
        %sum(current) + sum(abs(1-current))<=512
    %minimize (norm( resp_mat*current - resp_vec,2))
    minimize (-resp_vec'*log(RC+0.001) - (1-resp_vec)'*log(1-RC + 0.001))
cvx_end

resp_observed = resp_mat*current;


figure;
stem(resp_observed);
hold on;
stem(resp_vec);
legend('Observed','Desired');

figure;
stem(current);
xlabel('Electrode number');
ylabel('Current');


%% plot mosaic?


figure;
for icellType = 1:length(datarun.cell_types)
    subplot(3,5,icellType);
    plot_rf_fit(datarun, datarun.cell_types{icellType}.name,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
    for icell=1:length(final_selected_cells)
        if(sum(final_selected_cells(icell)==datarun.cell_types{icellType}.cell_ids)==1)
            plot_rf_fit(datarun,final_selected_cells(icell),'fill_color',[0,1,1],'alpha',0.3,'fill',true,'edge',true,'labels',false,'edge_color',[1,0,0]);
        end
    end
end

% colors for observed cells
% figure;
% for icell=1:nCells
% plot_rf_fit(datarun,cellIDs(icell),'fill_color',[0,1,1],'alpha',resp_observed(icell),'fill',true,'edge',true,'labels',false,'edge_color',[1,0,0]);
% end

responding_cells = cellIDs(resp_observed>0.01);
figure;
for icellType = 1:2%length(datarun.cell_types)
    subplot(2,1,icellType);
    plot_rf_fit(datarun, datarun.cell_types{icellType}.name,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
    for icell=1:length(responding_cells)
        if(sum(responding_cells(icell)==datarun.cell_types{icellType}.cell_ids)==1)
          
            plot_rf_fit(datarun,responding_cells(icell),'fill_color',[0,1,1],'alpha',resp_observed(cellIDs == responding_cells(icell)),'fill',true,'edge',true,'labels',false,'edge_color',[1,0,0]);
        end
    end
    
     for icell=1:length(final_selected_cells)
        if(sum(final_selected_cells(icell)==datarun.cell_types{icellType}.cell_ids)==1)
            plot_rf_fit(datarun,final_selected_cells(icell),'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[0,0,0],'sd_radius',1.2);
        end
     end
    axis equal
    xlim([10,70]);ylim([5,38]);
    set(gca,'xTick',[]);
    set(gca,'yTick',[]);
   title(datarun.cell_types{icellType}.name) 
end
