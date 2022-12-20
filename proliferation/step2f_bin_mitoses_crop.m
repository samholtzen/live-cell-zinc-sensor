%this one takes the signals and plots them aligned to mitosis to visualize
%average FRET curve after mitosis
rng(1)

for c=conditions_to_plot
    
    if ismember(c, [1,2])
        continue
    end
    
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
    
    CFP_store = CFP_store(:,1:180);
    YFP_store = YFP_store(:,1:180);
    H2B_store = H2B_store(:,1:180);
    for ii = 1:numel(mitosis_store)
        mitosis_store{ii} = mitosis_store{ii}(mitosis_store{ii} < 180);
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
            elseif i > 61 && i < 120
                FRET_align_second = [FRET_align_second;FRET_store_temp];
            elseif i > 120
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
    
    subplot(1,3,1)
    plot(frame_align./5, smooth(FRET_mean_first),'Color',colors_cell{c},'DisplayName', 'Early Mitoses','LineWidth',2)
    axis([-5 11 y_min y_max])
    resting_FRET_before = mean(FRET_mean_first(num_frames-20:num_frames-10));
    yline(resting_FRET_before,'--','Color',colors_cell{c},'LineWidth',2)
    xline(0,'--','DisplayName','Mitosis')
    title(['\fontsize{20}','Early Mitoses'])
    hold on
    
    subplot(1,3,2)
    plot(frame_align./5, smooth(FRET_mean_second),'Color',colors_cell{c},'DisplayName', 'Mid Mitoses','LineWidth',2)
    axis([-5 11 y_min y_max])
    resting_FRET_before = mean(FRET_mean_second(num_frames-20:num_frames-10));
    yline(resting_FRET_before,'--','Color',colors_cell{c},'LineWidth',2)
    xline(0,'--','DisplayName','Mitosis')
    title(['\fontsize{20}','Mid Mitoses'])
    hold on
    
    subplot(1,3,3)
    plot(frame_align./5, smooth(FRET_mean_third),'Color',colors_cell{c},'DisplayName', 'Late Mitoses','LineWidth',2)
    axis([-5 11 y_min y_max])
    resting_FRET_before = mean(FRET_mean_third(num_frames-20:num_frames-10));
    yline(resting_FRET_before,'--','Color',colors_cell{c},'LineWidth',2)
    xline(0,'--','DisplayName','Mitosis')
    title(['\fontsize{20}','Late Mitoses'])
    hold on

    
    
end

