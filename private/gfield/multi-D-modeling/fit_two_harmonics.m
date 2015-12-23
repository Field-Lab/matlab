function [fit_amplitudes, fit_error] = fit_two_harmonics(signal, time_points, fundamental_period, initial_amplitudes)

%subtract mean from data
%signal = signal - mean(signal);
%signal = signal ./ ext(signal);

f = @(amplitudes)fit_two_sign_waves(amplitudes, fundamental_period, time_points, signal);

[fit_amplitudes, fit_error] = fminsearch(f, initial_amplitudes);

if 1
    signal_fit = two_sine_waves(fit_amplitudes, fundamental_period, time_points);
    figure(1)
    plot(signal_fit, 'r', 'LineWidth', 3)
    hold on
    plot(signal, 'k')
    hold off
end



% fit error
function fit_error = fit_two_sign_waves(amplitudes, fundamental_period, time_points, signal)


fit_output = two_sine_waves(amplitudes, fundamental_period, time_points);
    
fit_error = sqrt(mean((fit_output - signal).^2));   

if 0
    signal_fit = two_sine_waves(amplitudes, fundamental_period, time_points);
    figure(1)
    plot(signal_fit, 'r')
    hold on
    plot(signal, 'k')
    hold off
    pause
end



    
    
% fit function
function signal_output = two_sine_waves(amplitudes, fundamental_period, time_points)

t_vals = 1:1:time_points;

signal_output = amplitudes(1) * sin((2*pi*t_vals./fundamental_period) + amplitudes(3))...
        + amplitudes(2) * sin((4*pi*t_vals./(fundamental_period)) + 1.57);

signal_output(signal_output < 0) = 0;
    
    

    
    
