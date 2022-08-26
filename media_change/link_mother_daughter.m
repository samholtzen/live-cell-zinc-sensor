function [all_full_traces_CFP, all_full_traces_YFP, all_full_traces_H2B, all_mitosis,num_frames] = link_mother_daughter(well,sensor)

%extract data one track at a time

% create structures to save data
num_frames = length(well{1}.ellipse_id);
all_full_traces_CFP = nan(0, num_frames); % each row represents a full trace
all_full_traces_YFP = nan(0, num_frames); % each row represents a full trace
all_full_traces_H2B = nan(0, num_frames); % each row represents a full trace

all_mitosis = cell(0, 1); % each element stores the mitosis frame ID.

% create genealogy: each element stores the mother Cell ID
genealogy = nan(length(well), 1);
for i=1:length(well)
    genealogy(cell2mat(well{i}.daughters)) = i;
end

% find full traces
% strategy: search backwards in time. Find cells that present at the last
% frame. Then for each such cell, trace to mother, grand mother, etc, until
% reaching the first frame.
for i=1:length(well)
    if isnan(well{i}.ellipse_id(end)) % not present at the last frame
        continue;
    end
    
    % store the information for the current full trace
    curr_id = i; curr_full_trace_CFP = nan(1, num_frames); curr_mitosis = []; curr_full_trace_YFP = nan(1, num_frames);curr_full_trace_H2B = nan(1, num_frames);
    
    % iterate
    while 1
        % cell is present between Frame first_id and Frame last_id
        first_id = find(~isnan(well{curr_id}.ellipse_id), 1, 'first');
        last_id = find(~isnan(well{curr_id}.ellipse_id), 1, 'last');
        % change the following line for other signals
        if strcmpi(sensor, 'NES')
            curr_full_trace_CFP(first_id:last_id) = smooth(well{curr_id}.CFP_cytoring_mean(first_id:last_id));
            curr_full_trace_YFP(first_id:last_id) = smooth(well{curr_id}.YFP_cytoring_mean(first_id:last_id));
        elseif strcmpi(sensor, 'NLS')
            curr_full_trace_CFP(first_id:last_id) = smooth(well{curr_id}.CFP_nuc_mean(first_id:last_id));
            curr_full_trace_YFP(first_id:last_id) = smooth(well{curr_id}.YFP_nuc_mean(first_id:last_id));
        end
        curr_full_trace_H2B(first_id:last_id) = smooth(well{curr_id}.H2B_nuc_mean(first_id:last_id));
        
        % add mitosis Frame ID
        if (last_id ~= num_frames)
            curr_mitosis = [last_id, curr_mitosis];
        end
        
        % if reaching Frame 1, record the trace, exit
        if (first_id == 1)
            all_full_traces_CFP = cat(1, all_full_traces_CFP, curr_full_trace_CFP);
            all_full_traces_YFP = cat(1, all_full_traces_YFP, curr_full_trace_YFP);
            all_full_traces_H2B = cat(1, all_full_traces_H2B, curr_full_trace_H2B);
            
            all_mitosis = cat(1, all_mitosis, {curr_mitosis});
            break;
        end
        
        % otherwise, search mother. Exit if there is no mother
        curr_id = genealogy(curr_id);
        if isnan(curr_id)
            break;
        end
    end
end
