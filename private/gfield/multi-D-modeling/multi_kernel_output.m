function output_sigs = multi_kernel_output(offsets, generator_matrix)

[num_kernels, num_frames] = size(generator_matrix);

% force the stc nonlinearities to be centered at a generator signal of zero
stc_kernel_x_offset = 0;

output_sigs = zeros(num_kernels, num_frames);


for dm = 1:num_kernels
    if dm == 1
        temp_sigs= generator_matrix(dm,:); % the generator signals associated with the sta
        %temp_output = (temp_sigs - offsets(1)).^2 + offsets(2);
        %temp_output(temp_sigs < offsets(1)) = offsets(2);
        
        temp_output = (temp_sigs - offsets(1)).^2;
        temp_output(temp_sigs < offsets(1)) = 0;

        
        output_sigs(dm,:) = temp_output;
        
        
        
%         offsets(1)
%         offsets(2)
%         plot(temp_sigs, temp_output, 'k.')
%         pause
    else
        %output_sigs(dm,:) = (generator_matrix(dm,:)-stc_kernel_x_offset).^2 + offsets(3);
        output_sigs(dm,:) = (generator_matrix(dm,:)-stc_kernel_x_offset).^2 + 0;

    end
end


