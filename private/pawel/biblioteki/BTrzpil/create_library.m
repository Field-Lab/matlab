function [pulse_library,pulse_coordinates]=create_library(array_amplitude_range,array_current_range,list_relatives_amplitude,array_option)


%array_amplitude_range - all amplitude pulse
%array_current_range - 
%list_relatives_amplitude
%array_option
%pulse_library - pulse matrix
%pulse_coordinates - pulse starting points in library
    pulse_library=create_pulse(array_amplitude_range(1), array_current_range(1), cell2mat(list_relatives_amplitude(1)), array_option(1));
    pulse_coordinates=1;
    
%library=horzcat(pulse1,pulse2);
%for i=2:N


l=length(array_amplitude_range);
for i=2:l
    pulse_i=create_pulse(array_amplitude_range(i), array_current_range(i), cell2mat(list_relatives_amplitude(i)), array_option(i));
    pulse_library=horzcat(pulse_library,pulse_i);
    [m,n]=size(cell2mat(list_relatives_amplitude(i-1)))
    k(i-1)=(n+2)*3 %previous pulse length
    pulse_coordinates=horzcat(pulse_coordinates,sum(k)+1)
end
[m, n]=size(pulse_library)
pulse_coordinates=horzcat(pulse_coordinates, n+1)
end
    