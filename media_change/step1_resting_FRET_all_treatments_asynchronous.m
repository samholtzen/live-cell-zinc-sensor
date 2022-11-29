%takes in combined signals and plots mean FRET in asynchronous cells
%(surrogate for resting FRET, this variable gets used in every script downstream)

for c=conditions_to_plot
    
    condition_index = find(conditions_to_plot == c);
    
    %loop through the columns (conditions)
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
    
    %iterate through the tracks in each condition
    num_frames = size(CFP_store,2);
    
    current_FRET = YFP_store./CFP_store;
    current_mitosis = mitosis_store;
    
    filter_vec_2 = max(current_FRET,[],2,'omitnan') < FRET_max & min(current_FRET,[],2,'omitnan') > FRET_min;
    
    FRET_filter = current_FRET(filter_vec_2,:);
    mitosis_filter = current_mitosis(filter_vec_2,:);
    
    mitosis_exist = ~cellfun(@isempty,mitosis_filter);
    
    filtered_FRET = FRET_filter(mitosis_exist,:);
    filtered_mitosis = mitosis_filter(mitosis_exist,:);
    
    %create frame vectors
    frame_vec = 1:num_frames;
    
    %take the average FRET ratio to plot
    
    FRET_mean = mean(filtered_FRET,1,'omitnan');
    
    if isempty(FRET_mean)
        continue
    end
    SEM = std(filtered_FRET,0,1,'omitnan')/sqrt(size(filtered_FRET,1));
    [FRET_maximum,ind_max] = max(FRET_mean(:,media_change:end),[],2,'omitnan');
    [FRET_minimum,ind_min]= min(FRET_mean(:,media_change:end),[],2,'omitnan');
    
    %     curve1 = FRET_mean + SEM;
    %     curve2 = FRET_mean - SEM;
    %     x2 = [frame_vec, fliplr(frame_vec)];
    %     inBetween = [curve1, fliplr(curve2)];
    %     fill((x2-media_change)/5, inBetween, 'g','FaceAlpha',0.2, 'FaceColor',colors_cell{c}, 'EdgeAlpha',0);
    %     hold on
    plot((frame_vec-media_change)./5, smooth(FRET_mean),'Color',colors_cell{c},'LineWidth',2,'DisplayName',condition_cell{c})
    axis([-5 (num_frames-media_change)/5 y_min y_max])
    hold on
    

end

xline(0,'--','DisplayName','Media Change')
xlabel('\fontsize{16}Time (hours)')
ylabel('\fontsize{16}Mean FRET Ratio')
title(['\fontsize{20}Resting FRET - ' upper(cell_type)])
legend('Location','bestoutside')
ax=gca;
ax.PositionConstraint = 'innerposition';
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
