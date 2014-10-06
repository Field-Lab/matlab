function fit_params = fit_static_NL(spike_times, gen_signals)


% FIT NONLINEARITY


params.fit = 'exp';

switch params.fit
    case 'exp'

%        % get spike times
%        % if there are N spikes per bin, list that bin N times
%         spike_times = [];
%         for nn = 1:full(max(spikes(:,cc)))
%             spike_times = [spike_times; find( spikes(:,cc)>(nn-1) )];
%         end

        % define function that gives likelihood
        L = @(x)identify_SNL_scalars_fn(x,spike_times,gen_signals);

        % turn gradient on
        options = optimset('GradObj','on');

        % identify best a and b
        fit = fminunc(L,[0 0],options)';
%        fit = fminunc(L,[0 0])';

        fit_params.a = fit(1);
        fit_params.b = fit(2);
        fit_params.type = params.fit;
end
    






function [L,g] = identify_SNL_scalars_fn(x,spike_times,gen_signal)

% likelihood
L =  length(spike_times)*x(2) + x(1)*sum(gen_signal.*spike_times) - ...
    exp(x(2))*sum(exp(x(1)*gen_signal))  ;
% gradient
if nargout > 1
    g(1) = sum(gen_signal.*spike_times) - exp(x(2))*sum(gen_signal.*exp(x(1)*gen_signal));
    g(2) = length(spike_times) - exp(x(2))*sum(exp(x(1)*gen_signal));
end

L = -L;
g = -g;

