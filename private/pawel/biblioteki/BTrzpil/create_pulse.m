function pulse=create_pulse(amplitude_total, current_range, relatives_amplitude, option)

    %amplitude_range - (for example 1uA)
    %current_range - operating range (for example: 1uA,4uA,...)
    %relatives_amplitude - vector relative aplitude pulse (for example [2 2-3 -3 1 1])
    %option - artifact reduction option on the channel stimulated
    state_before=[0;0;1;0];
    state_connect=[1,0,1,0];
    state_after=[0;0;1;0];
    relative_amplitude_max=max(abs(relatives_amplitude));
    length_vector=length(relatives_amplitude);
    pulse=zeros(4,3*(length_vector+2));
    pulse(:,1)=state_before;
    pulse(:,2:3)=dac(amplitude_total,current_range,relatives_amplitude(1),relative_amplitude_max);
    for i=1:length(relatives_amplitude)
       pulse(:,1+3*i)=state_connect;
       pulse(:,2+3*i:3+3*i)=dac(amplitude_total,current_range,relatives_amplitude(i),relative_amplitude_max);
    end
    pulse(:,length(pulse)-2)=state_after;
    pulse(:,length(pulse)-1:length(pulse))=dac(amplitude_total,current_range,relatives_amplitude(length(relatives_amplitude)),relative_amplitude_max);
    pulse;
end