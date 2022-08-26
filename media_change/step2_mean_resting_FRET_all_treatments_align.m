

for c=conditions_to_plot
    
    condition_index = find(conditions_to_plot == c);
    
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
    
    first_mitoses = [];
    first_mit_FRET = [];
    
    min_YFP = min(YFP_store,[],2,'omitnan');
    min_CFP = min(CFP_store,[],2,'omitnan');
    
    current_FRET = YFP_store./CFP_store;
    current_mitosis = mitosis_store;
    
    filter_vec_2 = max(current_FRET,[],2,'omitnan') < FRET_max & min(current_FRET,[],2,'omitnan') > FRET_min;
    
    FRET_filter = current_FRET(filter_vec_2,:);
    mitosis_filter = current_mitosis(filter_vec_2,:);
    
    mitosis_exist = ~cellfun(@isempty,mitosis_filter);
    
    filtered_FRET = FRET_filter(mitosis_exist,:);
    filtered_mitosis = mitosis_filter(mitosis_exist,:);
    
    frame_vec = 1:num_frames;
    FRET_align_before = [];
    FRET_align_after = [];
    for track = 1:length(filtered_mitosis)
        
        %get the signals out from the "store" variables
        mitoses = filtered_mitosis{track};
        FRET = filtered_FRET(track,:);
        
        
        %if it divides, we can align it to mitosis
        for i=mitoses
            
            first_mit = i;
            
            FRET_store_temp = nan(1,2*num_frames);
            
            if i < media_change

                FRET_store_temp(num_frames-first_mit:end-first_mit-1)=FRET;
                FRET_align_before = [FRET_align_before;FRET_store_temp];
                
            elseif i > media_change
                
                FRET_store_temp(num_frames-first_mit:end-first_mit-1)=FRET;
                FRET_align_after = [FRET_align_after;FRET_store_temp];

            end
        end
        
        
    end
    
    frame_align = ((1:2*num_frames) - num_frames);
    
    FRET_mean_before = mean(FRET_align_before,1,'omitnan');
    FRET_mean_after = mean(FRET_align_after,1,'omitnan');
    
    if isempty(FRET_mean_after) || isempty(FRET_mean_before)
        continue
    end
    
    subplot(2,round(length(conditions_to_plot)/2),condition_index)
    plot(frame_align./5, smooth(FRET_mean_before),'Color',colors_cell{3},'LineWidth',2, 'DisplayName', [condition_cell{3}])
    hold on
    plot(frame_align./5, smooth(FRET_mean_after),'Color',colors_cell{c},'LineWidth',2, 'DisplayName', [condition_cell{c}])
    axis([-5 15 y_min y_max])
    legend()
        
    
    xline(0,'--','DisplayName','Mitosis')
    
    title(['\fontsize{20}Resting FRET - ' upper(cell_type)])
    xlabel('\fontsize{12}Time (hours)')
    ylabel('\fontsize{16}Mean FRET Ratio')
    legend()
end


