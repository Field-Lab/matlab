% testClassifier
% calculate the confusion matricies of emile's classifier versus the manual
% classification 

clear all
close all
clc

[number_class, types_class] = textread('/Volumes/Lab/Projects/mining/code/cellTypeClassifier/2005-07-26-0data000/prediction.txt', '%d %s');

[number_real, types_real, junk] = textread('/Volumes/Analysis/2005-07-26-0/data000/rf-summary-classification.txt', '%d %s %s');

parasol_real = strfind(types_real,'parasol');
midget_real = strfind(types_real,'midget');
large_real = strfind(types_real,'large');
SBC_real = strfind(types_real,'SBC');
amacrine_real = strfind(types_real,'amacrine');
crap_real = strfind(types_real,'crap');
edge_real = strfind(types_real,'edge');
contam_real = strfind(types_real,'contaminated');
dup_real = strfind(types_real,'duplicates');


parasol_real_mat=[];
midget_real_mat=[];
large_real_mat=[];
SBC_real_mat=[];
amacrine_real_mat=[];
crap_real_mat = [];
contam_real_mat = [];

dup_real_mat = [];
edge_real_mat = [];

for i =1:size(number_real,1)
    if ~isempty(parasol_real{i})
        parasol_real_mat(i,1) = number_real(i);
    end
    if ~isempty(midget_real{i})
        midget_real_mat(i,1) = number_real(i);
    end
    if ~isempty(large_real{i})
        large_real_mat(i,1) = number_real(i);
    end
    
    if ~isempty(SBC_real{i})
        SBC_real_mat(i,1) = number_real(i);
    end
    
    if ~isempty(amacrine_real{i})
        amacrine_real_mat(i,1) = number_real(i);
    end
    if ~isempty(crap_real{i})
        crap_real_mat(i,1) = number_real(i);
    end
        if ~isempty(contam_real{i})
        contam_real_mat(i,1) = number_real(i);
        end
        if ~isempty(dup_real{i})
        dup_real_mat(i,1) = number_real(i);
        end
        if ~isempty(edge_real{i})
       edge_real_mat(i,1) = number_real(i);
    end
end
p_true = parasol_real_mat(parasol_real_mat~=0);
m_true = midget_real_mat(midget_real_mat~=0);
l_true = large_real_mat(large_real_mat~=0);
sbc_true = SBC_real_mat(SBC_real_mat~=0);
a_true = amacrine_real_mat(amacrine_real_mat~=0);
c_true = crap_real_mat(crap_real_mat~=0);
contam_true = contam_real_mat(contam_real_mat~=0);
edge_true = edge_real_mat(edge_real_mat~=0);
dup_true = dup_real_mat(dup_real_mat~=0);


%% Classifier

parasol_class = strfind(types_class,'Parasol');
midget_class = strfind(types_class,'Midget');
large_class = strfind(types_class,'Large');
SBC_class = strfind(types_class,'SBC');
amacrine_class = strfind(types_class,'Amacrine');
crap_class = strfind(types_class,'other');


parasol_class_mat=[];
midget_class_mat=[];
large_class_mat=[];
SBC_class_mat=[];
amacrine_class_mat=[];
crap_class_mat = [];

for i =1:size(number_class,1)
    if ~isempty(parasol_class{i})
        parasol_class_mat(i,1) = number_class(i);
    end
    if ~isempty(midget_class{i})
        midget_class_mat(i,1) = number_class(i);
    end
    if ~isempty(large_class{i})
        large_class_mat(i,1) = number_class(i);
    end
    
    if ~isempty(SBC_class{i})
        SBC_class_mat(i,1) = number_class(i);
    end
    
    if ~isempty(amacrine_class{i})
        amacrine_class_mat(i,1) = number_class(i);
    end
    if ~isempty(crap_class{i})
        crap_class_mat(i,1) = number_class(i);
    end
    
end
p_class = parasol_class_mat(parasol_class_mat~=0);
m_class = midget_class_mat(midget_class_mat~=0);
l_class = large_class_mat(large_class_mat~=0);
sbc_class = SBC_class_mat(SBC_class_mat~=0);
a_class = amacrine_class_mat(amacrine_class_mat~=0);
c_class = crap_class_mat(crap_class_mat~=0);




% order = amacrine, large, midget, parasol, SBC, crap
confusion = nan(6,6);
confusion(1,1) = sum(ismember(a_true, a_class));
confusion(2,1) = sum(ismember(l_true, a_class));
confusion(3,1) = sum(ismember(m_true, a_class));
confusion(4,1) = sum(ismember(p_true, a_class));
confusion(5,1) = sum(ismember(sbc_true, a_class));
confusion(6,1) = sum(ismember(c_true, a_class));

confusion(1,2) = sum(ismember(a_true, l_class));
confusion(2,2) = sum(ismember(l_true, l_class));
confusion(3,2) = sum(ismember(m_true, l_class));
confusion(4,2) = sum(ismember(p_true, l_class));
confusion(5,2) = sum(ismember(sbc_true, l_class));
confusion(6,2) = sum(ismember(c_true, l_class));

confusion(1,3) = sum(ismember(a_true, m_class));
confusion(2,3) = sum(ismember(l_true, m_class));
confusion(3,3) = sum(ismember(m_true, m_class));
confusion(4,3) = sum(ismember(p_true, m_class));
confusion(5,3) = sum(ismember(sbc_true, m_class));
confusion(6,3) = sum(ismember(c_true, m_class));

confusion(1,4) = sum(ismember(a_true, p_class));
confusion(2,4) = sum(ismember(l_true, p_class));
confusion(3,4) = sum(ismember(m_true, p_class));
confusion(4,4) = sum(ismember(p_true, p_class));
confusion(5,4) = sum(ismember(sbc_true, p_class));
confusion(6,4) = sum(ismember(c_true, p_class));

confusion(1,5) = sum(ismember(a_true, sbc_class));
confusion(2,5) = sum(ismember(l_true, sbc_class));
confusion(3,5) = sum(ismember(m_true, sbc_class));
confusion(4,5) = sum(ismember(p_true, sbc_class));
confusion(5,5) = sum(ismember(sbc_true, sbc_class));
confusion(6,5) = sum(ismember(c_true, sbc_class));

confusion(1,6) = sum(ismember(a_true, c_class));
confusion(2,6) = sum(ismember(l_true, c_class));
confusion(3,6) = sum(ismember(m_true, c_class));
confusion(4,6) = sum(ismember(p_true, c_class));
confusion(5,6) = sum(ismember(sbc_true, c_class));
confusion(6,6) = sum(ismember(c_true, c_class));

total = sum(confusion(:))
confusion = confusion./total*100;

fmt=[repmat('%8.2f ',1,6) '\n'];
fprintf(fmt,confusion');
