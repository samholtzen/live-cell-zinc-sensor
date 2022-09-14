%takes in combined signals and plots mean FRET in asynchronous cells
%(surrogate for resting FRET, this variable gets used in every script downstream)
resting_FRET_downstream = [];

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
    
    num_frames = size(CFP_store,2);
    
    %Calculate FRET
    current_FRET = YFP_store./CFP_store;
    
    %filter by appropriate FRET ratio and store the FRET for
    %math downstream
    
    filter_vec = max(current_FRET,[],2,'omitnan') < FRET_max & min(current_FRET,[],2,'omitnan') > FRET_min;
    
    FRET_filter = current_FRET(filter_vec,:);
    mitosis_filter = mitosis_store(filter_vec,:);
    
    mitosis_exist = ~cellfun(@isempty,mitosis_filter);
    
    filtered_FRET = FRET_filter(mitosis_exist,:);
    filtered_mitosis = mitosis_filter(mitosis_exist,:);
    
    %create frame vectors
    frame_vec = 1:num_frames;
    
    %take the average FRET ratio to plot
    FRET_mean = mean(filtered_FRET,1,'omitnan');
    
    if isempty(FRET_mean)
        %If there aren't enough tracks, skip it
        continue
    end
    
%     SEM = std(filtered_FRET,0,1,'omitnan')/sqrt(size(filtered_FRET,1));
%     
%     
%     curve1 = FRET_mean + SEM;
%     curve2 = FRET_mean - SEM;
%     x2 = [frame_vec, fliplr(frame_vec)];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2/5, inBetween, 'g','FaceAlpha',0.2, 'FaceColor',colors_cell{c}, 'EdgeAlpha',0);
%     hold on
    plot(frame_vec./5, smooth(FRET_mean),'Color',colors_cell{c},'LineWidth',2,'DisplayName',condition_cell{c})
    hold on
    
    resting_FRET_downstream = [resting_FRET_downstream mean(FRET_mean(60:120))];
    
end

%change axes range
axis([0 num_frames/5 y_min 7.5])
xlabel('\fontsize{16}Time (hours)')
ylabel('\fontsize{16}Mean FRET Ratio (YFP/CFP)')
title(['\fontsize{20}Resting FRET - ' upper(cell_type)])
legend('Location','bestoutside')
ax = gca;
ax.FontSize = 16;

