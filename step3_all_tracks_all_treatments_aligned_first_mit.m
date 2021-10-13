%almost identical to step3 script, except instead of plotting averages, we
%are plotting individual tracks

for c=1:column_number
    
    FRET_align = [];
    resting_FRET_condition = resting_FRET_downstream(c);
    FRET_mean_store = [];
    CFP_store = [];
    YFP_store = [];
    H2B_store = [];
    mitosis_store = [];
    
    num_conditions = length(conditions_to_plot);
    counter = 0;

    
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
        plotted_tracks = 0;
        
        mean_H2B_max = mean(max(H2B_store,[],2));
        
        
        for track = 1:length(mitosis_store)
            
            frame_vec = 1:num_frames;
            %make the timestamp vector
            
            mitoses = mitosis_store{track};
            
            current_H2B = H2B_store(track,:);
            FRET = YFP_store(track,:)./CFP_store(track,:);
            H2B_norm = current_H2B/max(current_H2B);
            
            if max(FRET) <= resting_FRET_condition+1 && min(FRET) >= resting_FRET_condition-1
                
                if ~isempty(mitoses)
                    
                    first_mit = mitoses(1);
                    
                    frame_align_hours = (frame_vec - first_mit)./5;
                    
                    if ismember(c,conditions_to_plot)
                        subplot(2,ceil(num_conditions/2),find(conditions_to_plot==c))
                        plot(frame_align_hours, FRET,'Color', [0.7 0.7 0.7])
                        hold on
                    end
                    
                end
                
            end
        end
        
    end
    
end

for c=1:column_number
    
    if ismember(c,conditions_to_plot)
        counter=counter+1;
        subplot(2,ceil(num_conditions/2),find(conditions_to_plot==c))
        hold on
        plot((-174:175)./5, all_FRET_mean_nes(c,:),'Color',colors_cell{c},'LineWidth',3)
        title({condition_cell{c}, ' FRET Ratio Aligned to Mitosis'})
        xlabel('Hours Relative to Mitosis')
        ylabel('FRET Ratio')
        axis([-4 17 3 7])
    end
end