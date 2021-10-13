%this one takes the combined signals variable signals_nes_dead and sums all
%conditions together

clearvars -except signals_nes signals_nes_dead resting_FRET_downstream resting_FRET_dead_downstream all_FRET_mean_nes all_FRET_mean_nes_dead

all_FRET_mean_nes = [];
for c=1:6
    FRET_align = [];
    resting_FRET_condition = resting_FRET_downstream(c);

    for nd = 1:size(signals_nes,1)
        
        wells = signals_nes(nd,c,1);
        
        CFP_store = [];
        YFP_store = [];
        H2B_store = [];
        mitosis_store = [];
        for j=1:numel(wells)
            
            %cycling through wells
            
            well = wells{j};
            
            
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
        FRET_align_well = [];

        for track = 1:length(mitosis_store)
            
            frame_vec = 1:num_frames;
            %make the timestamp vector
            
            mitoses = mitosis_store{track};
            current_H2B = H2B_store(track,:);
            FRET = YFP_store(track,:)./CFP_store(track,:);
            H2B_norm = current_H2B/max(current_H2B);
            
            if max(FRET) <= resting_FRET_condition+1 && min(FRET) >= resting_FRET_condition-1
                %filtering using a typical fret ratio based on the DR of the
                %sensor
                
                if isempty(mitoses)
                    
                    [TF, prom] = islocalmax(H2B_norm, 'MinSeparation',60,'SamplePoints',frame_vec,'MinProminence',0.9);
                    
                    TF_ind = find(TF);
                    
                    if ~isempty(TF_ind) && current_H2B(TF_ind) >= mean_H2B_max
                        
                        mitoses = TF_ind;
                        
                    end
                    
                end
                
                if ~isempty(mitoses)
                    %if it divides, we can align it to mitosis
                    
                    first_mit = 50;
                    
                    FRET_store_temp = zeros(1,2*num_frames);
                    %for math
                    frame_align = frame_vec - first_mit;
                    
                    FRET_store_temp(num_frames-first_mit:end-first_mit-1)=FRET;
                    FRET_align_well = [FRET_align_well;FRET_store_temp];
                    
                end
                
            end
        end
        
        FRET_align = [FRET_align;FRET_align_well];
    
    end
    
    FRET_nums = FRET_align>0;
    FRET_num = sum(FRET_nums);
    
    FRET_sum = sum(FRET_align);
    FRET_mean = smooth(FRET_sum ./ FRET_num);
    
    if c==1
        plot((-174:175)./5, FRET_mean, '-r','LineWidth',2,'DisplayName','ZD')
        hold on
    elseif c==2
        plot((-174:175)./5, FRET_mean,'Color',[0.9290 0.6940 0.1250],'LineWidth',2,'DisplayName','ZD + 2µM Zn')
        hold on
    elseif c==3
        plot((-174:175)./5, FRET_mean, '-c','LineWidth',2,'DisplayName','MM')
        hold on
    elseif c==4
        plot((-174:175)./5, FRET_mean, '-g','LineWidth',2,'DisplayName','ZR 30')
        hold on
    elseif c==5
        plot((-174:175)./5, FRET_mean, 'Color',[0.8500 0.3250 0.0980],'LineWidth',2,'DisplayName','ZR 50')
        hold on
    elseif c==6
        plot((-174:175)./5, FRET_mean,'Color',[0.5 0.5 0.5],'LineWidth',2,'DisplayName','Non Chelex Serum')
        hold on
    end
    
    all_FRET_mean_nes = [all_FRET_mean_nes; FRET_mean'];
end

axis([-10 10 3.5 6.5])

xline(0, '--','DisplayName','Mitosis');
legend