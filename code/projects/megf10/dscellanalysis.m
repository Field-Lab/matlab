function[ds_struct] = dscellanalysis(NumSpikesCell, StimComb)

%Calculates vector averages and DS index for each spatial and temporal DG period for each cell   

%Inputs: NumSpikesCell: Average total number of spikes fired at each trial type for all cells

         %StimComb: All the stimulus combinations (col 1: spatial period, col 2: temporal period, col 3: direction)
         
%Outputs: rho: N temporal periods x M spatial periods where each cell has normalized spike numbers in each direction 

          %RHO: N temporal periods x M spatial periods where each cell has spike numbers in each direction 
          
          %Theta: Angles for each cell at each spatial and temporal period
          
          %Num: List of all temporal periods
          
          %DSindex: DS indices calculated for each cell each S & T period
          
          %Mag: Normalized vector sum for each S & T period
          
          %Magmax: Maximum Mag over all temporal periods for each spatial period
          
          %Magave: Average Mag for each spatial period
          
          %Angle: Direction of Mag for each S & T period
          
  %Sneha Ravi 
  %Last revision: 12-18-2012
          
          

spperiods = unique(StimComb(:,1));
mag = cell(length(unique(StimComb(:,2))), length(unique(StimComb(:,1))));
MAG = cell(length(unique(StimComb(:,2))), length(unique(StimComb(:,1))));
dsindex = cell(length(unique(StimComb(:,2))), length(unique(StimComb(:,1))));
angle = cell(length(unique(StimComb(:,2))), length(unique(StimComb(:,1))));
rho = cell(length(unique(StimComb(:,2))), length(unique(StimComb(:,1))));
RHO = cell(length(unique(StimComb(:,2))), length(unique(StimComb(:,1))));
theta = cell(length(unique(StimComb(:,2))), length(unique(StimComb(:,1))));
U = cell(length(unique(StimComb(:,2))), length(unique(StimComb(:,1))));
V = cell(length(unique(StimComb(:,2))), length(unique(StimComb(:,1))));

for i = 1:length(spperiods)
    
    inn = ismember(StimComb(:,1), spperiods(i))';
    NumSpikesCell2 = NumSpikesCell(:,inn);
    StimComb2 = StimComb(inn',:);
    
[rho1 theta1 num] = get_rhotheta(NumSpikesCell2, StimComb2);
RHO1 = rho1;
%[rho1 theta1] = normalize_zerodirection(rho1, theta1);
[rho1 theta1] = normalize_maxspikes(rho1, theta1);
%[rho1 theta1] = normalize_totalspikes(rho1, theta1);

[X Y] = get_XY(rho1, theta1);
[XX YY] = get_XY(RHO1, theta1);

%[U1 V1 angle1 mag1] = vector_average(X, Y);     %Calculate vector average
[U1 V1 angle1 mag1] = vector_sum(X, Y); %Calculate vector sum
[~, ~, ~, MAG1] = vector_sum(XX, YY);
%[U1 V1 angle1 mag1] = vector_max(rho1, theta1);     %Calculate vector max

% figure(1)
% close all;
% [T R Unew Vnew] = polar_plots_one(rho1,theta1,U1,V1,num,cellnum);
% close;
% raster_plots(StimComb,datarun,num,T,R,Unew, Vnew, cellnum)

figure();
polar_plots_all(U1,V1,num, mag1);
%handle = title(spperiods(i));
%set(handle,'Position',[50,2]);

mag(:,i) = mag1;
angle(:,i) = angle1;
MAG(:,i) = MAG1;

magmax(i,:) = max_allspeeds(mag1);
magave(i,:) = max_avespeeds(mag1);

null = nullmin(rho1);
%null = nullopp(rho1, angle1, theta1);

%dsindex1 = [];
dsindex1 = DS_index_one(mag1, null);
%dsindex1 = DS_index_two(mag1, null);

dsindex(:,i) = dsindex1;

%figure(2)
%plot_hists(mag1, dsindex1, magmax(i,:), magave(i,:))
RHO(:,i) = RHO1; 
rho(:,i) = rho1;
theta(:,i) = theta1;
U(:,i)= U1;
V(:,i)= V1;
end
ds_struct = v2struct(mag, MAG, dsindex, magmax, magave, angle, rho, RHO, theta, num, U, V);
end