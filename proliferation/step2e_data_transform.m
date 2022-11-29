
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
    
    for track = 1:numel(filtered_mitosis)
        
        %get the signals out from the "store" variables
        mitoses = filtered_mitosis{track};
        FRET = filtered_FRET(track,:);
        
        
        %if it divides, we can align it to mitosis
        for i=mitoses
            
            first_mit = i;
            FRET_store_temp = nan(1,2*num_frames);
            
            %create a zero matrix that's double the number of
            %frames
            
            %picomolar for 5900 pM
            % dynamic range is 1.55-1.7, use this to get upper and lower
            % limits of range
            FRET_store_temp(num_frames-first_mit:end-first_mit-1)=FRET;
            
            if c < 3
                
                if i > 10
                    FRET_align = [FRET_align;FRET_store_temp];
                end
                
            else
                
                if i > 140
                    FRET_align = [FRET_align;FRET_store_temp];
                end
            end
            
        end
        
        
    end
    
    frame_align = ((1:2*num_frames) - num_frames);
    
    %store these for each well
    FRET_mean = mean(FRET_align,1,'omitnan');
    
    if c==1
        r_min = 4.13;
        r_max_mid = r_min*1.55;
        r_max_upper = r_min*1.45;
        r_max_lower = r_min*1.65;
    end
    
    
    if isempty(FRET_mean)
        continue
    end
    
    %Conversion to zinc
    zinc_mean = 5300 * ((FRET_mean((num_frames-20):(num_frames+75)) - r_min)/(r_max_mid-r_min)).^(1/0.29);
    
    resting_zinc_before = mean(zinc_mean(1:10));
    
    curve1 = 5300 * ((FRET_mean((num_frames-20):(num_frames+75)) - r_min)/(r_max_upper-r_min)).^(1/0.29);
    curve2 = 5300 * ((FRET_mean((num_frames-20):(num_frames+75)) - r_min)/(r_max_lower-r_min)).^(1/0.29);
    x2 = [frame_align((num_frames-20):(num_frames+75)), fliplr(frame_align((num_frames-20):(num_frames+75)))];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2/5, inBetween, 'g','FaceAlpha',0.2, 'FaceColor',colors_cell{c}, 'EdgeAlpha',0);
    hold on;
    
    plot(frame_align((num_frames-20):(num_frames+75))./5, smooth(zinc_mean),'Color',colors_cell{c},'DisplayName', condition_cell{c},'LineWidth',2)
    hold on
    yline(resting_zinc_before,'--','Color',colors_cell{c},'LineWidth',2)
    title([{'\fontsize{18}'}; {upper(cell_type)}])
    xlabel('\fontsize{12}Time (hours)')
    ylabel('\fontsize{16}[Zn^2^+] (pM)')
    ax = gca;
    ylim([0,3000])
    ax.FontSize = 16;
    
end



