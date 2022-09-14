%this script plots one individual trace of a zinc spike after mitosis

rng(1)
close all
clearvars -except signals_nes signals_nes_dead resting_FRET_downstream resting_FRET_dead_downstream

for c=1

    FRET_align = [];
    for nd = 1:size(signals_nes,1)

        wells = signals_nes(nd,c,1);

        CFP_store = [];
        YFP_store = [];
        H2B_store = [];
        mitosis_store = [];
        for j=1:numel(wells)

            %cycling through wells

            well = wells{j};


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
                    curr_full_trace_CFP(first_id:last_id) = smooth(well{curr_id}.CFP_cytoring_mean(first_id:last_id));
                    curr_full_trace_YFP(first_id:last_id) = smooth(well{curr_id}.YFP_cytoring_mean(first_id:last_id));
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
            CFP_store = [CFP_store;all_full_traces_CFP];
            YFP_store = [YFP_store;all_full_traces_YFP];
            H2B_store = [H2B_store;all_full_traces_CFP];
            mitosis_store = [mitosis_store; all_mitosis];
        end
        plotted_tracks = 0;

        mean_H2B_max = mean(max(H2B_store,[],2));

        random = randi(length(mitosis_store),1,1);

        FRET_mean_store = [];
        for track = 1:length(mitosis_store)

            frame_vec = 1:num_frames;
            %make the timestamp vector

            mitoses = mitosis_store{track};

            current_H2B = H2B_store(track,:);
            FRET = YFP_store(track,:)./CFP_store(track,:);
            H2B_norm = (current_H2B-min(current_H2B))/(max(current_H2B)-min(current_H2B));

            if max(FRET) <= 4.2+1 && min(FRET) >= 4.2-1

                if ~isempty(mitoses)

                    FRET_mean_store = [FRET_mean_store; FRET];

                    first_mit = mitoses(1);

                    frame_align_hours = (frame_vec - first_mit)./5;
                    if c==1 && track == 341
                        yyaxis left
                        plot(frame_align_hours, FRET)
                        yyaxis right
                        plot(frame_align_hours, H2B_norm,'k')
                        xlim([-7,5])
                        break
                    end

                end

            end
        end

    end


end
