%this one takes the signals and plots them aligned to mitosis to visualize
%average FRET curve after mitosis
rng(1)

for c=conditions_to_plot
    
    condition_index = find(conditions_to_plot == c);
    
    FRET_align_first = [];
    FRET_align_second = [];
    FRET_align_third = [];
    
    
    FRET_align = [];
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
    
    current_FRET = YFP_store./CFP_store;
    current_mitosis = mitosis_store;
    
    %filter by appropriate FRET ratio and store the FRET for
    %math downstream
    
    filter_vec = max(current_FRET,[],2,'omitnan') < FRET_max & min(current_FRET,[],2,'omitnan') > FRET_min;
    
    FRET_filter = current_FRET(filter_vec,:);
    mitosis_filter = current_mitosis(filter_vec,:);
    mitosis_exist = ~cellfun(@isempty,mitosis_filter);
    
    filtered_FRET = FRET_filter(mitosis_exist,:);
    filtered_mitosis = mitosis_filter(mitosis_exist,:);
    
    %create frame fector
    frame_vec = 1:num_frames;
    
    
    for track = 1:numel(filtered_mitosis)
        
        %get the signals out from the "store" variables
        mitoses = filtered_mitosis{track};
        FRET = filtered_FRET(track,:);
        
        
        %if it divides, we can align it to mitosis
        for i=mitoses
            
            first_mit = i;
            
            %create a zero matrix that's double the number of
            %frames
            
            
            FRET_store_temp = nan(1,2*num_frames);
            
            
            %replace the indices in the FRET_store_temp vector that
            %correspond to the FRET aligned to mitosis with the
            %FRET values
            FRET_store_temp(num_frames-first_mit:end-first_mit-1)=FRET;
            if i < 60
                FRET_align_first = [FRET_align_first;FRET_store_temp];
            elseif i > 80 && i < 160
                FRET_align_second = [FRET_align_second;FRET_store_temp];
            elseif i > 160
                FRET_align_third = [FRET_align_third;FRET_store_temp];
            end
            
            
        end
        
        
    end
    
    frame_align = ((1:2*num_frames) - num_frames);
    
    %store these for each well
    FRET_mean_first = mean(FRET_align_first,1,'omitnan');
    FRET_mean_second = mean(FRET_align_second,1,'omitnan');
    FRET_mean_third = mean(FRET_align_third,1,'omitnan');
    
    if isempty(FRET_mean_first) || isempty(FRET_mean_second) || isempty(FRET_mean_third)
        continue
    end
    
    subplot(2,round(length(conditions_to_plot)/2),condition_index)
    
    
    plot(frame_align./5, smooth(FRET_mean_first),'Color',[0.2,0.2,0.2],'DisplayName', 'Early Mitoses','LineWidth',2)
    hold on
    plot(frame_align./5, smooth(FRET_mean_second),'Color',[0.5,0.5,0.5],'DisplayName', 'Mid Mitoses','LineWidth',2)
    hold on
    plot(frame_align./5, smooth(FRET_mean_third),'Color',[0.7,0.7,0.7],'DisplayName', 'Late Mitoses','LineWidth',2)
    hold on
    axis([-5 15 y_min y_max])
    xline(0,'--','DisplayName','Mitosis')
    title(['\fontsize{20}',condition_cell{c}])
    legend()
    
    
end

