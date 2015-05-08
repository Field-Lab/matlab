%% Generate responses
% Calculate filter output for each sub-unit for each frame and calculate
% number of spikes for each frame-bin (binned response) .. So that would be
% used for STA calculation ? 



mov2=zeros(Filtdim1 ,Filtdim2,movieLen+Filtlen-1);
mov2(:,:,Filtlen:movieLen+Filtlen-1)=mov; % Append zeros before the movie
nTrials=1;
SubUnit_Response_test_movie_script

% Spike triggered sub-unit Input 
spikeTriggeredSubUnitInput(binnedResponses,cell_resp)