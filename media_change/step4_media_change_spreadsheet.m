
out_table = table();

for c=conditions_to_plot
    
    table_conditions = {};
    cell_FRET = [];
    CFP_store = [];
    YFP_store = [];
    mitosis_store = [];
    mitosis_timing = {};
    
    for nd = 1:size(struct_cell,1)
        curr_struct = struct_cell{nd,c,1};
        
        if ~isempty(curr_struct)
            if ~isempty(curr_struct.YFP)
                CFP_store = [CFP_store;curr_struct.CFP];
                YFP_store = [YFP_store;curr_struct.YFP];
                mitosis_store = [mitosis_store; curr_struct.mitosis];
            end
        end
    end
        
    current_FRET = YFP_store./CFP_store;
    current_mitosis = mitosis_store;
    
    filter_vec_2 = max(current_FRET,[],2,'omitnan') < FRET_max & min(current_FRET,[],2,'omitnan') > FRET_min;
    
    filtered_FRET = current_FRET(filter_vec_2,:);
    filtered_mitosis = current_mitosis(filter_vec_2,:);
    

    
    for i=1:length(filtered_mitosis)
        
        cell_FRET = [cell_FRET;mean(filtered_FRET(i,media_change-20:media_change-1))];
        mitosis_timing = [mitosis_timing;'before'];
        cell_FRET = [cell_FRET;mean(filtered_FRET(i,media_change:end))];
        table_conditions = [table_conditions;condition_cell{c}];
        mitosis_timing = [mitosis_timing;'after'];
        table_conditions = [table_conditions;condition_cell{c}];
        
    end
    
    all_zn_data = table(cell_FRET,mitosis_timing,table_conditions);
    out_table = [out_table;all_zn_data];
    
end
writetable(out_table, ['/Volumes/hard_drive_1/zinc_spike','/',experiment_date, '_zn_spike_media_change_', cell_type, '.csv'])
writetable(out_table, [output_str,'/',experiment_date, '_zn_spike_media_change_', cell_type, '.csv'])
