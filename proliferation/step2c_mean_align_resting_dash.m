
rng(1)

for c=conditions_to_plot
    
    
    
    condition_index = find(conditions_to_plot == c);
    
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
    
    if ifcrop
        % Cropping out the last 60 frames of the movie to remove sections
        % of tracks that were captured during contact inhibition
        CFP_store = CFP_store(:,1:180);
        YFP_store = YFP_store(:,1:180);
        H2B_store = H2B_store(:,1:180);
        
        for ii = 1:numel(mitosis_store)
            mitosis_store{ii} = mitosis_store{ii}(mitosis_store{ii} < 180);
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
    rand_tracks = randi(numel(filtered_mitosis),1,1000);
    
    
    for track = rand_tracks
        
        %get the signals out from the "store" variables
        mitoses = filtered_mitosis{track};
        FRET = filtered_FRET(track,:);
        
        
        %if it divides, we can align it to mitosis
        for i=mitoses
            
            first_mit = i;
            FRET_store_temp = nan(1,2*num_frames);
            
            %create a zero matrix that's double the number of
            %frames
            
            FRET_store_temp(num_frames-first_mit:end-first_mit-1)=FRET;
            
            if c < 3
                
                if i > 10
                    FRET_align = [FRET_align;FRET_store_temp];
                end
                
            else
                
                if i > 120
                    FRET_align = [FRET_align;FRET_store_temp];
                end
            end
            
        end
        
        
    end
    
    frame_align = ((1:2*num_frames) - num_frames);
    
    %store these for each well
    FRET_mean = mean(FRET_align,1,'omitnan');
    
    
    
    if isempty(FRET_mean)
        continue
    end
    
    resting_FRET_before = mean(FRET_mean(num_frames-20:num_frames-10));
    
    plot(frame_align./5, smooth(FRET_mean),'Color',colors_cell{c},'DisplayName', condition_cell{c},'LineWidth',2)
    hold on
    yline(resting_FRET_before,'--','Color',colors_cell{c},'LineWidth',2)
    axis([-4 11 4 6])
    title([{'\fontsize{18}Early Mitosis - '}; {upper(cell_type)}])
    xlabel('\fontsize{12}Time (hours)')
    ylabel('\fontsize{16}Mean FRET Ratio')
    legend('off')
    ax = gca;
    ax.FontSize = 16;
    
end



