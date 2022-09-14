%this one takes the signals and plots them aligned to mitosis to visualize
%average FRET curve after mitosis
all_FRET_mean = [];
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
    
    % Take the average of 500 of the tracks to save on time
    if length(filtered_mitosis) > 500
        
        rand_tracks = randi(length(filtered_mitosis),1,500);
        
    else
        
        rand_tracks = 1:length(filtered_mitosis);
        
    end
    
    for track = rand_tracks
        
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
            FRET_align = [FRET_align;FRET_store_temp];
            
        end
        
        
    end
    
    frame_align = ((1:2*num_frames) - num_frames);
    
    %store these for each well
    FRET_mean = mean(FRET_align,1,'omitnan');
    if isempty(FRET_mean)
        continue
    end
    
    plot(frame_align./5, smooth(FRET_mean),'Color',colors_cell{c},'DisplayName', condition_cell{c},'LineWidth',2)
    hold on
    axis([-5 15 y_min y_max])
    legend()
    
    all_FRET_mean  = [all_FRET_mean;FRET_mean];
    
end

%find the number of cells that exist at a given frame and sum them
%together to get n

xline(0,'--','DisplayName','Mitosis')
title(['\fontsize{20}Resting FRET - ' upper(cell_type)])
xlabel('\fontsize{12}Time (hours)')
ylabel('\fontsize{16}Mean FRET Ratio (YFP/CFP)')
legend('Location','bestoutside')
ax = gca;
ax.FontSize = 16;
