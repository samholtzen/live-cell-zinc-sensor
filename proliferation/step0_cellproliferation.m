
for c=conditions_to_plot
    
    condition_index = find(conditions_to_plot == c);
    
    proliferation_store = [];
    
    for nd = 1:size(struct_cell,1)
        curr_struct = struct_cell{nd,c,1};
        
        if ~isempty(curr_struct)
            if ~isempty(curr_struct.YFP)
                proliferation_store = [proliferation_store;curr_struct.proliferation];
                
            end
        end
    end
    
    num_frames = size(proliferation_store,2);
    frame_vec = 1:num_frames;
    proliferation_mean = mean(proliferation_store,1);
    
    % Calculate SEM
    SEM = std(proliferation_store,[],1)/sqrt(size(proliferation_store,1));
    
    % Plot the error bands first
    curve1 = proliferation_mean + SEM;
    curve2 = proliferation_mean - SEM;
    x2 = [frame_vec, fliplr(frame_vec)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2/5, inBetween, 'g','FaceAlpha',0.2, 'FaceColor',colors_cell{c}, 'EdgeAlpha',0);
    hold on;
    
    % Plot data
    plot(frame_vec./5, mean(proliferation_store./proliferation_store(:,1),1),'Color',colors_cell{c},'LineWidth',4,'DisplayName',condition_cell{c})
    hold on
    
    
end

%change axes range
xlim([0 num_frames/5])


% Add labels and title
xlabel('\fontsize{16}Time (hours)')
ylabel('\fontsize{16}Normalized Cell Count')
title(['\fontsize{20}Proliferation - ' upper(cell_type)])
legend('Location','bestoutside')
ax = gca;
ax.FontSize = 16; 
