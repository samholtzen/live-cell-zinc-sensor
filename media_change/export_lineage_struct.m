%% Read in Files
% Navigate to the file that was saved from the "combine_folders.m" script
% and select the "signals.mat" file it saved.

% [filename,path] = uigetfile();
% 
% cell_type = input('What cell type is this? Input a string.');
% date = input('What date was this experiment conducted? Input a number in the format YYYYMMDD.');
% 
% load(fullfile(path,filename))

%% Initialize the storage cell array

num_rows = size(all_signals,1);
num_cols = size(all_signals,2);
struct_cell = cell(size(all_signals));



for i = 1:num_rows
    
    for j= 1:num_cols
        
        if ~isempty(all_signals{i,j})
            
            temp_struct = struct();
            
            % Change the details below to match the names of your
            % fluorophores.
            
            [all_CFP, all_YFP, all_H2B, all_mitosis, num_frames, well_ids]...
                = link_two(all_signals(i,j));
            
            proliferation = get_proliferation(all_signals(i,j));
            

            
            temp_struct.YFP = all_YFP;
            temp_struct.CFP = all_CFP;
            temp_struct.H2B = all_H2B;
            temp_struct.mitosis = all_mitosis;
            temp_struct.proliferation = proliferation;
            
            struct_cell{i,j} = temp_struct;
            
        end
        
    end
    
end

% save(fullfile(path,[num2str(date),'_',cell_type,'_struct_cell.mat']),'struct_cell')