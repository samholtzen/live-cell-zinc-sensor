%takes in combined signals and plots mean FRET in asynchronous cells
%(surrogate for resting FRET, this variable gets used in every script downstream)
resting_FRET_downstream = [];

for c = 1:6
    %1 is ZD, 2 is ZD + 2µM zinc, 3 is MM, 4 is ZR30
    
    FRET_mean_store = [];
    CFP_store = [];
    YFP_store = [];
    H2B_store = [];
    mitosis_store = [];
    for nd = 1:size(signals,1)
        
        wells = signals(nd,c,1);
        
        
        for j=1:numel(wells)
            
            %cycling through wells
            
            well = wells{j};
            
            [all_full_traces_CFP, all_full_traces_YFP, all_full_traces_H2B, all_mitosis,num_frames] = link_mother_daughter(well,sensor);
            
            CFP_store = [CFP_store;all_full_traces_CFP];
            YFP_store = [YFP_store;all_full_traces_YFP];
            H2B_store = [H2B_store;all_full_traces_CFP];
            mitosis_store = [mitosis_store; all_mitosis];
        end
        
        %iterate through the tracks in each condition
        for track = 1:length(mitosis_store)
            
            %Calculate FRET
            current_YFP = YFP_store(track,:);
            current_CFP = CFP_store(track,:);
            current_FRET = current_YFP./current_CFP;
            
            
            if max(current_FRET) <= 6 && min(current_FRET) >= 3
                
                %filter by appropriate FRET ratio and store the FRET for
                %math downstream
                FRET_mean_store = [FRET_mean_store; current_FRET];
                
            end
        end
    end
    
    %create frame fector
    frame_vec = 1:num_frames;
    
    %take the average FRET ratio to plot
    FRET_mean = mean(FRET_mean_store);
    
    %these statements will plot each of the conditions. change the labels
    %according to the conditions in each column, and comment out the
    %"plot..." line if you don't want it plotted
    
    if ismember(c,conditions_to_plot)
        plot(frame_vec./5, FRET_mean,'Color',colors_cell{c},'LineWidth',2,'DisplayName',condition_cell{c})
        hold on
    end
    
    resting_FRET_downstream = [resting_FRET_downstream mean(FRET_mean(60:120))];
    
end

%change axes range
axis([0 35 3 6])


xlabel('Time (hours)')
ylabel('Mean FRET Ratio (YFP/CFP)')
title('Resting FRET')
legend(condition_cell{conditions_to_plot})
