function [output_LOGI] = subR_compute_LOGISTIC(MAX,RATE, Y_INT, input_LOGI)
    OFFSET   = log( (MAX/Y_INT) - 1  ) / RATE;     
    output_LOGI = (MAX ./ (1 + exp(-RATE * (input_LOGI- OFFSET) )));
end