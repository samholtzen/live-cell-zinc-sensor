
out_table = table();

for c=conditions_to_plot
    
    table_conditions = {};
    resting_FRET = [];
    peak_FRET = [];
    mitosis_FRET = [];
    mitosis_frame = [];
    CFP_store = [];
    YFP_store = [];
    cell_FRET = [];
    FRET_drop = [];
    peak_index = [];
    mitosis_store = [];
    
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
    num_frames = size(CFP_store, 2);
    
    first_mitoses = [];
    first_mit_FRET = [];

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
    
    for i=1:length(filtered_mitosis)
        
        mitoses = filtered_mitosis{i};
        
        if mitoses(1) > 20 && mitoses(end) < num_frames-20
            
            %if mitosis events happen multiple times, we are
            %cycling through them
            
            for k=mitoses
                
                FRET_after = diff(filtered_FRET(i,k:k+20));
                
                FRET_drop = [FRET_drop;min(FRET_after)]; 
                cell_FRET = [cell_FRET; mean(filtered_FRET(i,:),'omitnan')];
                resting_FRET = [resting_FRET; mean(filtered_FRET(i,k-20:k-10),'omitnan')];
                [max_FRET, max_ind] = max(filtered_FRET(i,k:k+10),[], 'omitnan','linear');
                peak_FRET = [peak_FRET; max_FRET];
                peak_index = [peak_index;(max_ind-1)/5];
                mitosis_FRET = [mitosis_FRET;filtered_FRET(i,k)];
                mitosis_frame = [mitosis_frame;k];
                table_conditions = [table_conditions;condition_cell{c}];
            end
            
        end
        
    end
    %calculates zinc spike from the peak and resting FRET
    zn_spike = peak_FRET-resting_FRET;
    
    all_zn_data = table(cell_FRET,resting_FRET,peak_FRET,FRET_drop,mitosis_frame,peak_index,mitosis_FRET,zn_spike,table_conditions);
    out_table = [out_table;all_zn_data];
    
end

writetable(out_table, ['/Volumes/hard_drive_1/zinc_spike','/',experiment_date, '_zn_spike_', cell_type, '.csv'])
writetable(out_table, [output_str,'/',experiment_date, '_zn_spike_', cell_type, '.csv'])
