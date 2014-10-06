function dispsig(signalMatrix);


% DISPSIG - Display Signals
%
% DISPSIG(signalMatrix) displays a set of signals in a plot matrix.
% Assumes that the signals are column vectors of signalMatrix.
%
%
% EXAMPLE
%       dispsig (sig)
%       title('Signals')
%
%
% This function is needed by FASTICA and FASTICAG

% 16.3.1998

numSignals = size (signalMatrix, 2);

for i = 1:numSignals,
  subplot (numSignals,1, i);
  plot (signalMatrix (:,i));
end
% Move the handle to the first subplot, so that the title may
% easily be added to the top.
subplot (numSignals,1, 1);
