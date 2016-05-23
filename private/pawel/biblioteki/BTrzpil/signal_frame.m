function [frame, item]=signal_frame(item_frame)
frame0=frame_zero();
frame=repmat(frame0,1, item_frame);
item=item_frame;