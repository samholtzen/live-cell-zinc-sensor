%this one takes the signals and plots them aligned to mitosis to visualize
%average FRET curve after mitosis

%clears all the variables except the ones stored from previous scripts

all_FRET_mean_nes = [];
for c=1:column_number
    FRET_align = [];
    resting_FRET_condition = resting_FRET_downstream(c);
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
        
        plotted_tracks = 0;
        mean_H2B_max = mean(max(H2B_store,[],2));
        FRET_align_well = [];
        
        %cycle through each track in the well
        for track = 1:length(mitosis_store)
            
            %make the timestamp vector
            frame_vec = 1:num_frames;
            
            %get the signals out from the "store" variables
            mitoses = mitosis_store{track};
            current_H2B = H2B_store(track,:);
            FRET = YFP_store(track,:)./CFP_store(track,:);
            H2B_norm = current_H2B/max(current_H2B);
            
            if max(FRET) <= resting_FRET_condition+1 && min(FRET) >= resting_FRET_condition-1
                %filtering using a typical fret ratio based on the DR of the
                %sensor
                
                if isempty(mitoses)
                    %sometimes elliptrack can miss a mitosis event, so
                    %we are using a manual algorithm to find mitosis events
                    
                    %if a mitosis event occurs, H2B-Halo signal will spike
                    %at that point, which we can find using the
                    %"islocalmax" function
                    [TF, prom] = islocalmax(H2B_norm, 'MinSeparation',60,'SamplePoints',frame_vec,'MinProminence',0.9);
                    
                    TF_ind = find(TF);
                    
                    if ~isempty(TF_ind) && current_H2B(TF_ind) >= mean_H2B_max
                        
                        %if we find one we missed, add the mitosis events to the "mitosis" variable
                        %to use later
                        mitoses = TF_ind;
                        
                    end
                    
                end
                
                if ~isempty(mitoses)
                    %if it divides, we can align it to mitosis
                    
                    first_mit = mitoses(1);
                    mit_after_change = mitoses(mitoses > media_frame);
                    
                    %create a zero matrix that's double the number of
                    %frames
                    FRET_store_temp = zeros(1,2*num_frames);
                    
                    %align the frame vector to mitosis by subtracting
                    frame_align = frame_vec - first_mit;
                    
                    %replace the indices in the FRET_store_temp vector that
                    %correspond to the FRET aligned to mitosis with the
                    %FRET values
                    FRET_store_temp(num_frames-first_mit:end-first_mit-1)=FRET;
                    
                    %store this for math later
                    FRET_align_well = [FRET_align_well;FRET_store_temp];
                    
                end
                
            end
        end
        
        %store these for each well
        FRET_align = [FRET_align;FRET_align_well];
        
    end
    
    %find the number of cells that exist at a given frame and sum them
    %together to get n
    FRET_nums = FRET_align>0;
    FRET_num = sum(FRET_nums);
    
    %sum the FRET signals together
    FRET_sum = sum(FRET_align);
    
    %divide summed FRET signals by n cells to get the mean
    FRET_mean = smooth(FRET_sum ./ FRET_num);
    
    if ismember(c,conditions_to_plot)
        plot((-174:175)./5, FRET_mean,'Color',colors_cell{c},'LineWidth',2)
        hold on
    end
    
    all_FRET_mean_nes = [all_FRET_mean_nes; FRET_mean'];
end

axis([-10 10 3 6])

xline(0, '--','DisplayName','Mitosis');
legend(condition_cell{conditions_to_plot})