% Standalone function that is used exactly as main_file.m to plot how the
% mitosis events are distributed across the movie

clearvars

[filename,path] = uigetfile('MultiSelect','on');

if ~iscell(filename)
    filename = {filename};
end
S = load([path,filename{1}]);
fieldname = fieldnames(S);
struct_cell = S.(fieldname{1});
split_filename = split(filename{1},'_');

%Date of experiment
experiment_date = split_filename{1};

%Sensor used
sensor = 'nes';

cell_type = split_filename{2};

colors_cell = {'#33637a','#3DA5D9','#0C0808','#CE9600','#ED8A07','#B35C05'};

switch cell_type
    
    case {'scr','all'}
        condition_cell = {'ZD3','ZD2','Minimal Media','ZR15','ZR30','ZR50'};
        conditions_to_plot = [1,2,3,4,5,6];
        y_min = 4;
        y_max = 6.5;
        
    case 'kd'
        condition_cell = {'ZD3','ZD2','Minimal Media','ZR15','ZR30','ZR50'};
        conditions_to_plot = [1, 2, 3, 4, 5,6];
        y_min = 4;
        y_max = 6.5;
        
    case {'dead','wt'}
        colors_cell = {'#3DA5D9','#3DA5D9','#0C0808','#CE9600','#ED8A07','#B35C05'};
        condition_cell = {'ZD2','ZD + 3µM Zn','Minimal Media','ZR15','ZR30','No Chx'};
        conditions_to_plot = [1,2, 3, 4, 5,6];
        y_min = 3;
        y_max = 6.5;
        
end


output_str = [sensor,'_',cell_type,'_','output'];
if ~exist(output_str,'dir')
    mkdir(output_str)
end


% Run parameters for filtering YFP signal. These are empirically tested to
% ensure we remove cells that aren't expressing the construct.

FRET_min = 3;
FRET_max = 8;

for c=conditions_to_plot
    
    condition_index = find(conditions_to_plot == c);
    
    FRET_filter = [];
    CFP_store = [];
    YFP_store = [];
    H2B_store = [];
    mitosis_store = [];
    
    for nd = 1:size(struct_cell,1)
        curr_struct = struct_cell{nd,c,1};
        
        if ~isempty(curr_struct)
            if ~isempty(curr_struct.YFP)
                CFP_store = [CFP_store;curr_struct.CFP];
                YFP_store = [YFP_store;curr_struct.YFP];
                H2B_store = [H2B_store;curr_struct.H2B];
                mitosis_store = [mitosis_store; curr_struct.mitosis];
            end
        end
    end
    
    num_frames = size(CFP_store, 2);
    
    %initialize storage vectors
    first_mitoses = [];
    first_mit_FRET = [];
    
    current_FRET = YFP_store./CFP_store;
    current_mitosis = mitosis_store;
    
    %filter by appropriate FRET ratio and store the FRET for
    %math downstream
    filter_vec_2 = max(current_FRET,[],2,'omitnan') < FRET_max & min(current_FRET,[],2,'omitnan') > FRET_min;
    
    FRET_filter = current_FRET(filter_vec_2,:);
    mitosis_filter = current_mitosis(filter_vec_2,:);
    
    mitosis_diff = cellfun(@diff, mitosis_filter, 'UniformOutput',false);
    mitosis_diff_log = cellfun(@(x) x<40, mitosis_diff, 'UniformOutput', false);
    mitosis_diff_sum = cellfun(@sum, mitosis_diff_log) < 1;
    
    mitosis_exist = ~cellfun(@isempty,mitosis_filter);
    
    filtered_FRET = FRET_filter(mitosis_exist ,:);
    filtered_mitosis = mitosis_filter(mitosis_exist ,:);
    
    %create frame vector
    frame_vec = 1:num_frames;
    
    mit_arr = [];    
    
    for i=1:length(filtered_mitosis)
        
        mit_arr = [mit_arr, filtered_mitosis{i}(1:end-1)];
        
    end
    
    subplot(6,1,condition_index)
    h = histogram(mit_arr,50);
    h.FaceColor = colors_cell{c};
    h.EdgeColor = colors_cell{c};
    xlim([0,num_frames])
    xline(120,'--')
    ylim([0,100])
    hold on
    
end