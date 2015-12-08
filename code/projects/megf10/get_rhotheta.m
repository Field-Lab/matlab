function[rho theta num] = get_rhotheta(NumSpikesCell2, StimComb2)

%Function extracts spike numbers and angles of each direction for each temporal period

%Inputs: NumSpikesCell2: Total avergae spike numbers for each trial 

        %StimComb2: List of all combinations of S & T periods & directions
        
%Outputs: rho: spike numbers of each direction for each temporal period

          %theta: angles of each direction for each temporal period
          
          %num: List of all temporal periods
          
          %Sneha Ravi 
  %Last revision: 12-18-2012

num = unique(StimComb2(:,2)); %Number of temporal periods
rho = cell(length(num), 1);
theta= cell(length(num), 1);
for i = 1:length(num)
    [r,t] = deal([]);
    in = ismember(StimComb2(:,2), num(i))';
    r = NumSpikesCell2(:,in);
    t = (StimComb2(in,3)')*pi/180;
    t = repmat(t,size(r, 1), 1);
    rho{i,1} = r;
    theta{i,1} = t;
end