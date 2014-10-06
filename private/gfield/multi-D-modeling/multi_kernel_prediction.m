function firing_rate = multi_kernel_prediction(weights, offsets, generator_matrix)

output_sigs = multi_kernel_output(offsets, generator_matrix);

firing_rate = weights * output_sigs;



