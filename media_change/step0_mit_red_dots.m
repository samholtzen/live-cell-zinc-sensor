%This plots mitosis events as red dots over traces of FRET ratio in dark
%gray
rng(1)

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
    mitosis_diff_log = cellfun(@(x) x<0, mitosis_diff, 'UniformOutput', false);
    mitosis_diff_sum = cellfun(@sum, mitosis_diff_log) < 1;
    
    mitosis_exist = ~cellfun(@isempty,mitosis_filter);

    filtered_FRET = FRET_filter(mitosis_exist & mitosis_diff_sum,:);
    filtered_mitosis = mitosis_filter(mitosis_exist & mitosis_diff_sum,:);

    %create frame fector
    frame_vec = 1:num_frames;
    
    mit_arr = [];
    mit_FRET_arr = [];
    rand_tracks = randi(length(filtered_mitosis),1,100);
    
    for i=rand_tracks
        
        filtered_mit_before = filtered_mitosis{i}(1:end-1);
        
        mit_arr = [mit_arr, filtered_mit_before];
        mit_FRET_arr = [mit_FRET_arr, filtered_FRET(i,filtered_mit_before)];
        
    end
    
    if isempty(filtered_FRET)
        continue
    end
    
    % Plot all random track FRET ratios
    subplot(2,round(length(conditions_to_plot)/2),condition_index)
    plot(frame_vec./5, filtered_FRET(rand_tracks,:), 'Color', [0.25 0.25 0.25]);
    hold on
    
    %plot the stored mitosis events
    subplot(2,round(length(conditions_to_plot)/2),condition_index)
    plot(mit_arr./5,mit_FRET_arr, '.r','MarkerSize',10);
    title(['\fontsize{20}',condition_cell{c}])
    xlabel('\fontsize{12}Time (h)')
    ylabel('\fontsize{16}Mean FRET Ratio')
    xline(media_change/5,'--')
    axis([0 num_frames/5 FRET_min FRET_max])
    hold on
    
end





