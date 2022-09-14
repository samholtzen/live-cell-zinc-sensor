%takes in combined signals and plots mean FRET in asynchronous cells
%(surrogate for resting FRET, this variable gets used in every script downstream)
rng(1)
out_table = table();
for c=1:6
    coef_temp = [];
    condition_index = find(conditions_to_plot == c);
    
    %loop through the columns (conditions)
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
    
    %iterate through the tracks in each condition
    num_frames = size(CFP_store,2);
    
    %Calculate FRET
    
    
    %loop through the cells in one condition
    min_YFP = mean(YFP_store,2,'omitnan');
    min_CFP = mean(CFP_store,2,'omitnan');
    
    
    filter_vec_1 = min_CFP > 75;
    
    
    current_FRET = YFP_store(filter_vec_1,:)./CFP_store(filter_vec_1,:);
    current_mitosis = mitosis_store(filter_vec_1);
    %filter by appropriate FRET ratio and store the FRET for
    %math downstream
    
    filter_vec_2 = max(current_FRET,[],2,'omitnan') < FRET_max & min(current_FRET,[],2,'omitnan') > FRET_min;
    rand_tracks = randi(length(filtered_mitosis),1,500);
    
    filtered_FRET = current_FRET(filter_vec_2,:);
    mitosis_filter = current_mitosis(filter_vec_2,:);
    frame_vec = 1:num_frames;
    
    nan_find = isnan(filtered_FRET);
    nan_sum = sum(nan_find,1);
    
    if isempty(FRET_mean)
        continue
    end
    
    plot(frame_vec, nan_sum, '.')
    hold on

    
end