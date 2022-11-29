
rng(1)

out_cell = [];

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
    
    around_mitosis = (num_frames-20):(num_frames+75);
    
    %Conversion to zinc
    zinc_mean = 5300 * ((FRET_mean(around_mitosis) - r_min)/(r_max_mid-r_min)).^(1/0.29);
    zinc_max = 5300 * ((FRET_mean(around_mitosis) - r_min)/(r_max_upper-r_min)).^(1/0.29);
    zinc_min = 5300 * ((FRET_mean(around_mitosis) - r_min)/(r_max_lower-r_min)).^(1/0.29);
    
    resting_zinc_mean = mean(zinc_mean(1:10));
    resting_zinc_max = mean(zinc_max(1:10));
    resting_zinc_min = mean(zinc_min(1:10));
    
    peak_zinc_mean = max(zinc_mean(20:40),[], 'omitnan');
    peak_zinc_max = max(zinc_max(20:40),[], 'omitnan');
    peak_zinc_min = max(zinc_min(20:40),[], 'omitnan');
    
    resting_store = [resting_zinc_min; resting_zinc_mean; resting_zinc_max];
    peak_store = [peak_zinc_min; peak_zinc_mean; peak_zinc_max];
    
    out_cell = [out_cell;[resting_store,peak_store]];
    
    
end



