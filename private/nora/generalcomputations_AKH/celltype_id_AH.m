% UGLY CODE to find the index of datarun.cell_types .. in a fairly robust
% manner !! !!

% you get the idea %
% AK Heitman    start: 2013-08-28

% Inputs must be one of the following:
%        On-Parasol, Off-Parasol, On-Midget, Off-Midget, SBC
% 

% To minimize hassle of future coding

function [celltype_index, CIDs ] = celltype_id_AH(string_celltype, cell_dataruncelltypes)

ctypes = cell_dataruncelltypes; types = length(cell_dataruncelltypes);

match = zeros(1,types);  % this will turn to 1 when we have a match 
if strcmp(string_celltype , 'On Parasol')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'On Parasol') || strcmp(name, 'On Par') || strcmp(name, 'on parasol') ||...
                strcmp(name, 'on par') ||  strcmp(name, 'ON Parasol') ||  strcmp(name, 'On parasol')
            match(i_type) = 1;
            break
        end
    end
end

if strcmp(string_celltype , 'Off Parasol')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'Off Parasol') || strcmp(name, 'Off Par') || strcmp(name, 'off parasol') ||...
                strcmp(name, 'off par') ||  strcmp(name, 'OFF Parasol')||  strcmp(name, 'Off parasol')
            match(i_type) = 1;
            break
        end
    end
end

if strcmp(string_celltype , 'On Midget')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'On Midget') ||  strcmp(name, 'on midget') ||  strcmp(name, 'ON Midget')
            match(i_type) = 1;
            break
        end
    end
end

if strcmp(string_celltype , 'Off Midget')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'Off Midget') ||  strcmp(name, 'off midget') ||  strcmp(name, 'OFF Midget')
            match(i_type) = 1;
            break
        end
    end
end

if strcmp(string_celltype , 'SBC')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'SBC') ||  strcmp(name, 'sbc') ||  strcmp(name, 'Small Bistratified') ...
                ||strcmp(name, 'small bistratified')  
            match(i_type) = 1;
            break
        end
    end
end

if strcmp(string_celltype , 'On Large')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'ON Large') ||  strcmp(name, 'On Large') ...
             
            match(i_type) = 1;
            break
        end
    end
end
if strcmp(string_celltype , 'Off Large')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'OFF Large') ||  strcmp(name, 'Off Large') 
            match(i_type) = 1;
            break
        end
    end
end
if strcmp(string_celltype , 'crap')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'crap') ||  strcmp(name, 'Crap')   
            match(i_type) = 1;
            break
        end
    end
end
if strcmp(string_celltype , 'Unclassified')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'unclassified') ||  strcmp(name, 'Unclassified') 
            match(i_type) = 1;
            break
        end
    end
end

if strcmp(string_celltype , 'Crap')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'Crap') 
            match(i_type) = 1;
            break
        end
    end
end

if strcmp(string_celltype , 'nc5')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'nc5') 
            match(i_type) = 1;
            break
        end
    end
end

if strcmp(string_celltype , 'nc6')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'nc6') 
            match(i_type) = 1;
            break
        end
    end
end

if strcmp(string_celltype , 'nc4')
    
    for i_type = 1:types
        name = ctypes{i_type}.name;
        
        if strcmp(name, 'nc4') 
            match(i_type) = 1;
            break
        end
    end
end






celltype_index = find(match);
try
    CIDs           = ctypes{celltype_index}.cell_ids;
catch
    CIDs = [];
end

        
        
        
        
        
    
    
