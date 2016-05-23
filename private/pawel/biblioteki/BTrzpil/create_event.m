function event=create_event(pulse_library, pulse_coordinates, event_time, channel,time, pulse_number)

%pulse_library - pulse matrix
%pulse_coordinates - pulse starting points in library (matrix 1D)
%event_time - event diuration time
event_length=event_time/0.00005
event=repmat(frame_zero(),1,event_length);
%first_pulse_frame - pulse starting points in event

        first_pulse_frame=time/0.00005;
        
        pulse=pulse_library(:,pulse_coordinates(pulse_number):pulse_coordinates(pulse_number+1)-1);

        for i=1:length(pulse)/3
            pulse_block_start=(first_pulse_frame)*250+250*(i-1)+((channel*3)-2);
            pulse_block_end=(first_pulse_frame)*250+250*(i-1)+((channel*3));
            event(:,pulse_block_start:pulse_block_end)=pulse(:,(i*3-2):i*3);
        end
    

end