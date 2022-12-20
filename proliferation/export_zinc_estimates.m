
rng(1)
peak_store = [];
resting_store = [];
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
    
    rand_tracks = randi(numel(filtered_mitosis),1,2000);
    
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
                
                if i > 120
                    FRET_align = [FRET_align;FRET_store_temp];
                end
            end
            
        end
        
        
    end
    
    frame_align = ((1:2*num_frames) - num_frames);
    
    %store these for each well
    FRET_mean = mean(FRET_align,1,'omitnan');
    
    if c==1
        r_min = min(FRET_mean)-0.01;
        r_max_mid = r_min*1.55;
        r_max_upper = r_min*1.45;
        r_max_lower = r_min*1.65;
    end
    
    
    if isempty(FRET_mean)
        continue
    end
    
    FRET_to_convert = FRET_mean((num_frames-20):(num_frames+55));

    %Conversion to zinc
    zinc_mean = 5300 * ((FRET_to_convert - r_min)./(r_max_mid-FRET_to_convert)).^(1/0.29);
    zinc_mean_upper = 5300 * ((FRET_to_convert - r_min)./(r_max_upper-FRET_to_convert)).^(1/0.29);
    zinc_mean_lower = 5300 * ((FRET_to_convert - r_min)./(r_max_lower-FRET_to_convert)).^(1/0.29);
    
    before_index = 1:9;
    after_index = 21:32;
    
    resting = [mean(zinc_mean_lower(before_index));mean(zinc_mean(before_index));mean(zinc_mean_upper(before_index))];
    peak = [max(zinc_mean_lower(after_index));max(zinc_mean(after_index));max(zinc_mean_upper(after_index))];
    
    peak_store = [peak_store;peak];
    resting_store = [resting_store;resting];
    
end

spike = peak_store - resting_store;

out = [resting_store, peak_store, spike];
