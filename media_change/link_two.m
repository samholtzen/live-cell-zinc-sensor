function [all_FP1, all_FP2, all_NM, all_mitosis, num_frames, well_ids] = link_two(signals)

all_FP1 = [];
all_FP2 = [];
all_NM = [];
all_mitosis = [];
well_ids = [];

for c=1:numel(signals)
    
    wells = signals(c);
    
    
    for j=1:numel(wells)
        
        %cycling through wells
        
        cell_track = wells{j};
        
        
        %extract data one track at a time
        
        % create structures to save data
        num_frames = length(cell_track{1}.ellipse_id);
        all_full_traces_FP1 = nan(0, num_frames); % each row represents a full trace
        all_full_traces_FP2 = nan(0, num_frames); % each row represents a full trace
        all_full_traces_NM = nan(0, num_frames); % each row represents a full trace
        
        all_mitosis_temp = cell(0, 1); % each element stores the mitosis frame ID.
        
        % create genealogy: each element stores the mother Cell ID
        genealogy = nan(length(cell_track), 1);
        for i=1:length(cell_track)
            genealogy(cell2mat(cell_track{i}.daughters)) = i;
        end
        
        % find full traces
        % strategy: search backwards in time. Find cells that present at the last
        % frame. Then for each such cell, trace to mother, grand mother, etc, until
        % reaching the first frame.
        for i=1:length(cell_track)
            if isnan(genealogy(i)) % if the cell divides
                continue;
            end
            % store the information for the current full trace
            curr_id = i; curr_full_trace_FP1 = nan(1, num_frames); curr_mitosis = []; curr_full_trace_FP2 = nan(1, num_frames);curr_full_trace_NM = nan(1, num_frames);
            
            % iterate
            while 1
                % cell is present between Frame first_id and Frame last_id
                first_id = find(~isnan(cell_track{curr_id}.ellipse_id), 1, 'first');
                last_id = find(~isnan(cell_track{curr_id}.ellipse_id), 1, 'last');
                % change the following line for other signals
                curr_full_trace_FP1(first_id:last_id) = smooth(cell_track{curr_id}.CFP_cytoring_mean(first_id:last_id));
                curr_full_trace_FP2(first_id:last_id) = smooth(cell_track{curr_id}.YFP_cytoring_mean(first_id:last_id));
                curr_full_trace_NM(first_id:last_id) = smooth(cell_track{curr_id}.H2B_nuc_mean(first_id:last_id));
                
                % add mitosis Frame ID
                if (last_id ~= num_frames)
                    curr_mitosis = [last_id, curr_mitosis];
                end
                
                % if reaching Frame 1, record the trace, exit

                
                % otherwise, search mother. Exit if there is no mother
                curr_id = genealogy(curr_id);
                if isnan(curr_id)
                    all_full_traces_FP1 = cat(1, all_full_traces_FP1, curr_full_trace_FP1);
                    all_full_traces_FP2 = cat(1, all_full_traces_FP2, curr_full_trace_FP2);
                    all_full_traces_NM = cat(1, all_full_traces_NM, curr_full_trace_NM);
                    well_ids = [well_ids;c];
                    
                    all_mitosis_temp = cat(1, all_mitosis_temp, {curr_mitosis});
                    break;
                end
            end
        end
        all_FP1 = [all_FP1;all_full_traces_FP1];
        all_FP2 = [all_FP2;all_full_traces_FP2];
        all_NM = [all_NM;all_full_traces_FP1];
        all_mitosis = [all_mitosis; all_mitosis_temp];
    end
end

end