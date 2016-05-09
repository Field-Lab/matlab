function prediction = model_pred(b, reg_gen_signal_center, reg_gen_signal_surround)
    center = conv(reg_gen_signal_center, time_filter(b(1:6)), 'valid'); 
    surround = conv(reg_gen_signal_surround, time_filter(b(7:12)), 'valid'); 
    prediction = b(13)/(b(14)+(center+surround).^b(15))+b(16);
end