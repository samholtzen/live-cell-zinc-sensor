
rng(1)
out_table = table();
table_conditions = {};

for c=1:6
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
    
    num_frames = size(CFP_store,2);
    
    min_YFP = mean(YFP_store,2,'omitnan');
    min_CFP = mean(CFP_store,2,'omitnan');
    
    
    
    current_FRET = YFP_store./CFP_store;
    current_mitosis = mitosis_store;
    
    filter_vec_2 = max(current_FRET,[],2,'omitnan') < FRET_max & min(current_FRET,[],2,'omitnan') > FRET_min;
    
    filtered_FRET = current_FRET(filter_vec_2,:);
    
    if c==1
        r_min = min(mean(filtered_FRET,1,'omitnan'))-0.01;
        r_max_mid = r_min*1.55;
    end
    
    for i=1:20
        
        if length(filtered_FRET(:,1)) < 500
            rand_tracks = randi(length(filtered_FRET(:,1)),1,100);
        else
            rand_tracks = randi(length(filtered_FRET(:,1)),1,300);
        end
        
        FRET = filtered_FRET(rand_tracks,:);
        
        FRET_mean = mean(FRET,1,'omitnan');
        
        
        [FRET_maximum,ind_max] = max(FRET_mean(:,media_change:end),[],2,'omitnan');
        [FRET_minimum,ind_min]= min(FRET_mean(:,media_change:end),[],2,'omitnan');
        
        out_max = ind_max/5;
        out_min = ind_min/5;
        
        norm_FRET = FRET_mean(:,media_change:end) ./ FRET_mean(:,media_change);
        auc_FRET = trapz(norm_FRET,2);
        
        settle_FRET = mean(FRET_mean(end-30:end));
        
        mean_FRET_before = mean(FRET_mean(media_change-30:media_change-1),'all','omitnan');
        
        if settle_FRET >= r_max_mid || settle_FRET < r_min
            zinc_settle = NaN;
        else
            zinc_settle = 5300 * ((settle_FRET - r_min)./(r_max_mid-settle_FRET)).^(1/0.29);

        end
        
        zinc_before = 5300 * ((mean_FRET_before - r_min)./(r_max_mid-mean_FRET_before)).^(1/0.29);
        
        table_conditions = [table_conditions;condition_cell{c}];
        table_temp = table(FRET_maximum,FRET_minimum,out_max,settle_FRET,out_min,auc_FRET,zinc_before,zinc_settle);
        table_temp.table_conditions = repmat(condition_cell(c), size(table_temp.FRET_maximum));
        if ismember(c,[1,2,3])
            
            table_temp.t_half = NaN(size(table_temp.FRET_maximum));
            
        elseif ismember(c,[4,5,6])
            
            FRET_half_max = (mean(FRET_maximum,'omitnan') + mean_FRET_before)/2;
            
            mean_FRET_after = FRET_mean(media_change:end-30);
            
            FRET_after_sub = mean_FRET_after - FRET_half_max;
            
            FRET_pos = FRET_after_sub(ind_max:end) > 0;
            t_half = find(FRET_pos,1,'last')/5;
            
            table_temp.t_half = repmat(t_half,size(table_temp.FRET_maximum));
            
        end
        
        out_table = [out_table;table_temp];
        
    end
    
    
end

writetable(out_table, ['/Volumes/hard_drive_1/zinc_spike','/',experiment_date, '_additional_media_change_', cell_type, '.csv'])
writetable(out_table, [output_str,'/',experiment_date, '_additional_media_change_', cell_type, '.csv'])

