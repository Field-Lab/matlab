function sig_str = pop_motion_signal_colleen(velocity, spikes, indices1, indices2, x_pos, trigger, trial_length, tau, tol)
%POP_MOTION_SIGNAL calculates the cumulative motion signal for a given velocity between all pairs in population of neurons.
%   sig_str = POP_MOTION_SIGNAL(velocity, spikes, indices1, indices2, x_pos, trigger, stop, tau)
%   caluclates the signal strength for the given velocity where:
%
%       spikes -> is a [num_neurons x 1] cell containing [num_spikes x 1] 
%                 vectors corresponding to all the spike times of that
%                 specific cell.  Contains all the neurons picked up by the
%                 array.
%       indices1 -> contains the indices (for x_pos) for the neurons
%                   identified using white noise to check receptive fields
%       indices2 -> contains the indices (for spikes) to identify the
%                   neurons during each trial run.  indices1(i) and
%                   indices2(i) correspond to the same neuron.
%       x_pos -> is a [num_neurons2 x 1] vector containing the x position
%                of all the neurons identified in the white noise stimuli.
%                i.e. the x coordinate of the gaussian fit of the receptive
%                field.
%       trigger ->the start time for the trial you are using
%       trial_length -> length of the trial
%       tau -> the tuning parameter to use.  The width of the gaussians in
%              the filtered response algorithm
%
%   sig_str = POP_MOTION_SIGNAL(..., tol)
%   same as above with:
%       tol -> tolerance to use for the integration algorithm in matlab.
%              Defaults to 1e-3
%
%   For examples on usage see motion_script.m
%
% Marvin Thielk 2013
% mthielk@salk.edu
%
% Modified for speed by Malcolm Campbell 2014
% malcolmc@stanford.edu
% 
% Modified for even more speed by Colleen Rhoades 2014
% rhoades@stanford.edu

if nargin < 9
    tol = 1e-3;
end

sig_str = 0;

% create a list of all possible pairs between the neurons
  RelTol = tol;
    AbsTol = tol * .1;
    
       t = linspace(0, trial_length, 100);
%      valueAtEachT = zeros(length(t), length(indices2),2);
          valueAtEachT = zeros(length(t), 1,2);

     

for i = 1:1%length(indices2) % for every cell shift every other cell relative to it
        dx = x_pos(indices1) - x_pos(indices1(i));


  
    spks_1 = spikes{indices2(i)};
    spks_2 = spikes(indices2); % Cell array of spike trains of all types considering
    
    % trial start @ t=0
    spks_1 = spks_1 - trigger;
   
        spks_2=cellfun(@(x) x-trigger, spks_2,'UniformOutput',false);
%     spks_2{c} = spks_2{c} - trigger;
    
    % only consider spikes that occured in the trial
    spks_1 = spks_1(spks_1 >= 0 & spks_1 <= trial_length);
    spks_2 = cellfun(@(y) y(y >= 0 & y <= trial_length), spks_2, 'UniformOutput', false);

    
%     if abs(dx / velocity) > trial_length / 2 ||isempty(spks_1) || isempty(spks_2)
%     
%         flt_rsp1 = @(t) 0;
%         flt_rsp2 = @(t) 0;
%         flt_rsp1_shiftedRight = @(t) 0;
%         flt_rsp2_shiftedRight = @(t) 0;
%         flt_rsp1_shiftedLeft = @(t) 0;
%         flt_rsp2_shiftedLeft = @(t) 0;
%         spks_1_shiftedRight = [];
%         spks_2_shiftedRight = [];
%         spks_1_shiftedLeft = [];
%         spks_2_shiftedLeft = [];
%     return
%     end

    
    % circularly shift spikes by dt
    spks_1_shiftedRight = spks_1 - dx(i) / velocity;
    spks_1_shiftedRight = sort(mod(spks_1_shiftedRight, trial_length));
    
    spks_1_shiftedLeft = spks_1 + dx(i) / velocity;
    spks_1_shiftedLeft = sort(mod(spks_1_shiftedLeft, trial_length));
    
    dx_cell = num2cell(dx);
    spks_2_shiftedRight = cellfun(@(z, c) z - c / velocity, spks_2,dx_cell, 'UniformOutput', false); %right shift
    spks_2_shiftedRight = cellfun(@(z) sort(mod(z, trial_length)), spks_2_shiftedRight, 'UniformOutput', false);
    spks_2_shiftedLeft = cellfun(@(z, c) z + c / velocity, spks_2,dx_cell, 'UniformOutput', false); %right shift
    spks_2_shiftedLeft = cellfun(@(z) sort(mod(z, trial_length)), spks_2_shiftedLeft, 'UniformOutput', false);
    
    % replicate spikes before and after trial to minimize artifacts of spikes
    % shifting circularly across the border
    idx=cellfun('isempty',spks_2_shiftedLeft);
    spks_2_shiftedLeft(idx)={0}; %It replaces all empty cells with number 0
      
    idx=cellfun('isempty',spks_2_shiftedRight);
    spks_2_shiftedRight(idx)={0}; %It replaces all empty cells with number 0
    
    idx=isempty(spks_1_shiftedLeft);
    spks_1_shiftedLeft(idx)=0; %It replaces all empty cells with number 0idx=cellfun('isempty',spks_2_shiftedLeft);
       
    idx=isempty(spks_1_shiftedRight);
    spks_1_shiftedRight(idx)=0; %It replaces all empty cells with number 0idx=cellfun('isempty',spks_2_shiftedLeft);
       
    spks_1_shiftedRight = [spks_1_shiftedRight(ceil(end/2):end) - trial_length; spks_1_shiftedRight; spks_1_shiftedRight(1:floor(end/2)) + trial_length];
    spks_1_shiftedLeft = [spks_1_shiftedLeft(ceil(end/2):end) - trial_length; spks_1_shiftedLeft; spks_1_shiftedLeft(1:floor(end/2)) + trial_length];
    


    spks_2_shiftedRight = cellfun(@(z) [z(ceil(end/2):end) - trial_length; z; z(1:floor(end/2)) + trial_length], spks_2_shiftedRight, 'UniformOutput', false);
    spks_2_shiftedLeft = cellfun(@(z) [z(ceil(end/2):end) - trial_length; z; z(1:floor(end/2)) + trial_length], spks_2_shiftedLeft, 'UniformOutput', false);
    
    % filter responses
    flt_rsp1 = filtered_response(spks_1, tau);
    flt_rsp2 = cellfun(@(x) filtered_response(x, tau), spks_2, 'UniformOutput', false);
    flt_rsp1_shiftedRight = filtered_response(spks_1_shiftedRight, tau);
    flt_rsp2_shiftedRight = cellfun(@(x) filtered_response(x, tau), spks_2_shiftedRight, 'UniformOutput', false);
      flt_rsp1_shiftedLeft = filtered_response(spks_1_shiftedLeft, tau);
    flt_rsp2_shiftedLeft = cellfun(@(x) filtered_response(x, tau), spks_2_shiftedLeft, 'UniformOutput', false);
     
   
    
   
 
    for c= 1:length(t)
        time = t(c);
        % First column is every cell aligned to cell 1
        % First row is value of first pairing at t = 0
        [right, left] = summedCellsAtT(time, flt_rsp1, flt_rsp2, flt_rsp1_shiftedRight, flt_rsp2_shiftedRight, flt_rsp1_shiftedLeft, flt_rsp2_shiftedLeft, indices1);
        valueAtEachT(c,i,1) = right;
        valueAtEachT(c,i,2) = left;
    end
    
    
    


    
    
    
%         str = integral(@(t) summedCellsAtT(t, flt_rsp1, flt_rsp2, flt_rsp1_shifted, flt_rsp2_shifted, indices1)^2, 0, trial_length, ...
%         'AbsTol', AbsTol, 'RelTol', RelTol);
    %str = 0;
    
    
end

sig_str = sum(valueAtEachT(:,:,1).^2)-sum(valueAtEachT(2,:,:).^2);

end

