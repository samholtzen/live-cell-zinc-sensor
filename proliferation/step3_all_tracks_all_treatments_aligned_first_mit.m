%almost identical to step3 script, except instead of plotting averages, we
%are plotting individual tracks

%This plots mitosis events as red dots over traces of FRET ratio in dark
%gray
rng(1)
%columns are different media conditions, rows are replicates

for c=conditions_to_plot
    %loop through the columns (conditions)
    
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
    CFP_store = CFP_store(:,1:180);
    YFP_store = YFP_store(:,1:180);
    H2B_store = H2B_store(:,1:180);
    for ii = 1:numel(mitosis_store)
        mitosis_store{ii} = mitosis_store{ii}(mitosis_store{ii} < 180);
    end
    num_frames = size(CFP_store, 2);
    %% Find and join mother/daughter pairs to create full lineage
    
    
    %cycling through wells
    %% Get mitosis events and use those to plot FRET ratio and red dots
    
    %create random tracks to plot
    
    %initialize storage vectors
    first_mitoses = [];
    first_mit_FRET = [];
    
    %loop through the cells in one condition
    min_YFP = mean(YFP_store,2,'omitnan');
    min_CFP = mean(CFP_store,2,'omitnan');
    
    
    
    
    
    current_FRET = YFP_store./CFP_store;
    current_mitosis = mitosis_store;
    %filter by appropriate FRET ratio and store the FRET for
    %math downstream
    
    filter_vec_2 = max(current_FRET,[],2,'omitnan') < FRET_max & min(current_FRET,[],2,'omitnan') > FRET_min;
    
    FRET_filter = current_FRET(filter_vec_2,:);
    mitosis_filter = current_mitosis(filter_vec_2,:);
    
    mitosis_exist = ~cellfun(@isempty,mitosis_filter);
    
    
    
    
    filtered_FRET = FRET_filter(mitosis_exist,:);
    filtered_mitosis = mitosis_filter(mitosis_exist,:);
    
    
    %create frame fector
    frame_vec = 1:num_frames;
    if isempty(filtered_FRET)
        continue
    end
    
    rand_tracks = randi(length(filtered_mitosis),1,100);
    
    for i=rand_tracks
        
        first_mit = filtered_mitosis{i}(ceil(end/2));
        frame_mit = frame_vec - first_mit;
        
        %take the average FRET ratio to plot
        subplot(2,round(length(conditions_to_plot)/2),condition_index)
        plot(frame_mit./5, filtered_FRET(i,:), 'Color', [0.75 0.75 0.75]);
        hold on
        
    end
    
    subplot(2,round(length(conditions_to_plot)/2),condition_index)
    plot(((1:num_frames*2)-num_frames)./5, smooth(all_FRET_mean(condition_index,:)), 'Color',colors_cell{c},'LineWidth',5);
    axis([-5 11 FRET_min FRET_max])
    title(['\fontsize{20}',condition_cell{c}])
    xlabel('\fontsize{12}Time (h)')
    ylabel('\fontsize{16}Mean FRET Ratio')
    hold on
    
end

