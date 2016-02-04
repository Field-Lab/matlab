% [each_cell_type, each_acf] = get_just_timecourse(3);

for i = [1,3, 4, 5,6,7,15,19,20,21]
    if i == 1
        test_tc1 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{9}];
        
        
        test_acf1 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{9}]';
        label1 = [ones(size(each_cell_type{i}{1},2),1); 2*ones(size(each_cell_type{i}{3},2),1); 3*ones(size(each_cell_type{i}{9},2),1)];
        
    elseif i ==2
        test_tc2 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{9}, each_cell_type{i}{12}, each_cell_type{i}{10}, each_cell_type{i}{13}, each_cell_type{i}{8}];
        
        test_acf2 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{9}; each_acf{i}{12}; each_acf{i}{10}; each_acf{i}{13}; each_acf{i}{8}]';
        
        
        label2 = [21*ones(size(each_cell_type{i}{1},2),1); 22*ones(size(each_cell_type{i}{3},2),1); 23*ones(size(each_cell_type{i}{9},2),1); 24*ones(size(each_cell_type{i}{12},2),1); 25*ones(size(each_cell_type{i}{10},2),1); 26*ones(size(each_cell_type{i}{13},2),1); 27*ones(size(each_cell_type{i}{8},2),1)];
    elseif i ==3
        test_tc3 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{11}];
        
        
        test_acf3 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{11}]';
        label3 = [31*ones(size(each_cell_type{i}{1},2),1); 32*ones(size(each_cell_type{i}{3},2),1); 34*ones(size(each_cell_type{i}{11},2),1)];
        
    elseif i ==4
        test_tc4 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{11}, each_cell_type{i}{9}];
        
        for t = 1:size(test_tc4,2)
            test_tc4_temp(:,t) = interp1(1:30, test_tc4(:,t), 1:0.5:30)';
        end
        
        test_tc4 = test_tc4_temp(size(test_tc4_temp,1)-29:size(test_tc4_temp,1),:);
        test_acf4 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{11}; each_acf{i}{9}]';
        label4 = [41*ones(size(each_cell_type{i}{1},2),1); 42*ones(size(each_cell_type{i}{3},2),1); 43*ones(size(each_cell_type{i}{11},2),1); 44*ones(size(each_cell_type{i}{9},2),1)];
        
        
    elseif i ==5
        test_tc5 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{9}];
        
        %         for t = 1:size(test_tc5,2)
        %             test_tc4_temp(:,t) = interp1(1:30, test_tc4(:,t), 1:0.5:30)';
        %         end
        %
        %         test_tc4 = test_tc4_temp(size(test_tc4_temp,1)-29:size(test_tc4_temp,1),:);
        test_acf5 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{9}]';
        label5 = [51*ones(size(each_cell_type{i}{1},2),1); 52*ones(size(each_cell_type{i}{3},2),1); 53*ones(size(each_cell_type{i}{9},2),1)];
 
        
            elseif i ==6
        test_tc6 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{13}];
        
        %         for t = 1:size(test_tc5,2)
        %             test_tc4_temp(:,t) = interp1(1:30, test_tc4(:,t), 1:0.5:30)';
        %         end
        %
        %         test_tc4 = test_tc4_temp(size(test_tc4_temp,1)-29:size(test_tc4_temp,1),:);
        test_acf6 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{13}]';
        label6 = [61*ones(size(each_cell_type{i}{1},2),1); 62*ones(size(each_cell_type{i}{3},2),1); 63*ones(size(each_cell_type{i}{13},2),1)];
    
        
        
     elseif i ==7
        test_tc7 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{11}];
        
        %         for t = 1:size(test_tc5,2)
        %             test_tc4_temp(:,t) = interp1(1:30, test_tc4(:,t), 1:0.5:30)';
        %         end
        %
        %         test_tc4 = test_tc4_temp(size(test_tc4_temp,1)-29:size(test_tc4_temp,1),:);
        test_acf7 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{11}]';
        label7 = [71*ones(size(each_cell_type{i}{1},2),1); 72*ones(size(each_cell_type{i}{3},2),1); 73*ones(size(each_cell_type{i}{11},2),1)];
 
        
            
     elseif i ==15
        test_tc15 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{12}];
         for t = 1:size(test_tc15,2)
            test_tc15_temp(:,t) = interp1(1:30, test_tc15(:,t), 1:0.5:30)';
        end
        
        test_tc15 = test_tc15_temp(size(test_tc15_temp,1)-29:size(test_tc15_temp,1),:);
        %         for t = 1:size(test_tc5,2)
        %             test_tc4_temp(:,t) = interp1(1:30, test_tc4(:,t), 1:0.5:30)';
        %         end
        %
        %         test_tc4 = test_tc4_temp(size(test_tc4_temp,1)-29:size(test_tc4_temp,1),:);
        test_acf15 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{12}]';
        label15 = [151*ones(size(each_cell_type{i}{1},2),1); 152*ones(size(each_cell_type{i}{3},2),1); 153*ones(size(each_cell_type{i}{12},2),1)];
            
     elseif i ==19
        test_tc19 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{15}];
         for t = 1:size(test_tc19,2)
            test_tc19_temp(:,t) = interp1(1:30, test_tc19(:,t), 1:0.5:30)';
        end
        
        test_tc19 = test_tc19_temp(size(test_tc19_temp,1)-29:size(test_tc19_temp,1),:);
        %         for t = 1:size(test_tc5,2)
        %             test_tc4_temp(:,t) = interp1(1:30, test_tc4(:,t), 1:0.5:30)';
        %         end
        %
        %         test_tc4 = test_tc4_temp(size(test_tc4_temp,1)-29:size(test_tc4_temp,1),:);
        test_acf19 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{15}]';
        label19 = [191*ones(size(each_cell_type{i}{1},2),1); 192*ones(size(each_cell_type{i}{3},2),1); 193*ones(size(each_cell_type{i}{15},2),1)];
   
     elseif i ==20
        test_tc20 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{9}];
         for t = 1:size(test_tc20,2)
            test_tc20_temp(:,t) = interp1(1:30, test_tc20(:,t), 1:0.5:30)';
        end
        
        test_tc20 = test_tc20_temp(size(test_tc20_temp,1)-29:size(test_tc20_temp,1),:);
        %         for t = 1:size(test_tc5,2)
        %             test_tc4_temp(:,t) = interp1(1:30, test_tc4(:,t), 1:0.5:30)';
        %         end
        %
        %         test_tc4 = test_tc4_temp(size(test_tc4_temp,1)-29:size(test_tc4_temp,1),:);
        test_acf20 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{9}]';
        label20 = [201*ones(size(each_cell_type{i}{1},2),1); 202*ones(size(each_cell_type{i}{3},2),1); 203*ones(size(each_cell_type{i}{9},2),1)];
   
       
         elseif i ==21
        test_tc21 = [each_cell_type{i}{1}, each_cell_type{i}{3}, each_cell_type{i}{8}];
         for t = 1:size(test_tc21,2)
            test_tc21_temp(:,t) = interp1(1:30, test_tc21(:,t), 1:0.5:30)';
        end
        
        test_tc21 = test_tc21_temp(size(test_tc21_temp,1)-29:size(test_tc21_temp,1),:);
        %         for t = 1:size(test_tc5,2)
        %             test_tc4_temp(:,t) = interp1(1:30, test_tc4(:,t), 1:0.5:30)';
        %         end
        %
        %         test_tc4 = test_tc4_temp(size(test_tc4_temp,1)-29:size(test_tc4_temp,1),:);
        test_acf21 = [each_acf{i}{1}; each_acf{i}{3}; each_acf{i}{8}]';
        label21 = [211*ones(size(each_cell_type{i}{1},2),1); 212*ones(size(each_cell_type{i}{3},2),1); 213*ones(size(each_cell_type{i}{8},2),1)];
   
       
           
    
        
    end
    
    
end


test1 = [test_acf1; test_tc1];
% test2 = [test_acf2; test_tc2];
test3 = [test_acf3; test_tc3];
test4 = [test_acf4; test_tc4];
test5 = [test_acf5; test_tc5];
test6 = [test_acf6; test_tc6];
test7 = [test_acf7; test_tc7];
test15 = [test_acf15; test_tc15];
test19 = [test_acf19; test_tc19];
test20 = [test_acf20; test_tc20];
test21 = [test_acf21; test_tc21];

test = [test1, test3, test4, test5, test6, test7, test15, test19, test20, test21];

label = [label1;label3; label4;label5; label6; label7;label15; label19; label20;label21];



pairs = [ 1 2 3;
    31 32 33;
    41 42 43;
    51 52 53;
    61 62 63;
    71 72 73 ;
    151 152 153;
    191 192 193 ;
    201 202 203; 
    211 212 213];
    

for k = 1:size(pairs,1)
% on midget
divide_by = repmat(mean(test(:,label == pairs(k,1))'), size(test(:,label== pairs(k,2)), 2),1)';
divide_by(abs(divide_by) < 0.01) = 0.01;
test(:,label==pairs(k,2)) = test(:,label==pairs(k,2))./divide_by;
test(:,label==pairs(k,2)) = filter([1/3 1/3 1/3], 1, test(:,label==pairs(k,2)));



% on smooth
divide_by = repmat(mean(test(:,label == pairs(k,1))'), size(test(:,label==pairs(k,3)), 2),1)';
divide_by(abs(divide_by) < 0.01) = 0.01;
test(:,label==pairs(k,3)) = test(:,label==pairs(k,3))./divide_by;
test(:,label==pairs(k,3)) = filter([1/3 1/3 1/3], 1, test(:,label==pairs(k,3)));


end

% % on midget
% divide_by = repmat(mean(test(:,label == 41)'), size(test(:,label==42), 2),1)';
% divide_by(abs(divide_by) < 0.01) = 0.01;
% test(:,label==42) = test(:,label==42)./divide_by;
% test(:,label==42) = filter([1/3 1/3 1/3], 1, test(:,label==42));
% 
% % on smooth
% 
% divide_by = repmat(mean(test(:,label == 41)'), size(test(:,label==43), 2),1)';
% divide_by(abs(divide_by) < 0.01) = 0.01;
% test(:,label==43) = test(:,label==43)./divide_by;
% test(:,label==43) = filter([1/3 1/3 1/3], 1, test(:,label==43));
% % on large 2
% divide_by = repmat(mean(test(:,label == 41)'), size(test(:,label==44), 2),1)';
% divide_by(abs(divide_by) < 0.01) = 0.01;
% test(:,label==44) = test(:,label==44)./divide_by;
% test(:,label==44) = filter([1/3 1/3 1/3], 1, test(:,label==44));
% 
% 
% % on smooth
% 
% divide_by = repmat(mean(test(:,label ==51)'), size(test(:,label==53), 2),1)';
% divide_by(abs(divide_by) < 0.01) = 0.01;
% test(:,label==53) = test(:,label==53)./divide_by;
% test(:,label==53) = filter([1/3 1/3 1/3], 1, test(:,label==53));


% label = [label1;label2];

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(test(:,label==2 | label==3 | label ==32| label ==34| label ==42| label ==43| label ==52| label ==53| label ==62| label ==63| label ==72| label ==73| label ==152| label ==153| label ==192| label ==193| label ==202| label ==203|  label ==212 | label ==213)', 'NumComponents',3);

[COEFF, SCORE_tc, LATENT, TSQUARED, EXPLAINED, MU] = pca(test_tc', 'NumComponents',3);

[COEFF, SCORE_acf, LATENT, TSQUARED, EXPLAINED, MU] = pca(test_acf', 'NumComponents',3);

labels = label(label==2 | label==3 | label ==32| label ==34| label ==42| label ==43| label ==52| label ==53| label ==62| label ==63| label ==72| label ==73| label ==152| label ==153| label ==192| label ==193| label ==202| label ==203|  label ==212 | label ==213);
for i = 1:length(labels)
    a = num2str(labels(i));
    if strcmp(a(end), '2')
        labels_real(i) = 1;
    elseif strcmp(a(end), '3') 
        labels_real(i) = 2;
    elseif strcmp(a(end), '4')
         labels_real(i) = 3;

    end
end


figure; gscatter(SCORE(:,1), SCORE(:,2), label(label==2 | label==3 | label ==32| label ==34| label ==42| label ==43| label ==52| label ==53| label ==62| label ==63| label ==72| label ==73| label ==152| label ==153| label ==192| label ==193| label ==202| label ==203|  label ==212 | label ==213))


figure; gscatter(SCORE_tc(:,1), SCORE_tc(:,2), label)


