%this script plots individual tracks for NES sensor in each condition in a 2x3 subplot
%and uses the variable all_signals for one row

%added the resting fret from asynchronous to this, need to add it to dead
%too


for c=1:6
    
    resting_FRET_condition = resting_FRET_downstream(c);
    FRET_align = [];
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
        
        resting_FRET = [];
        peak_FRET = [];        
        %cycle through the cells in each well
        for track = 1:length(mitosis_store)
            
            frame_vec = 1:num_frames;
            %make the timestamp vector
            
            %get out mitosis event information
            mitoses = mitosis_store{track};
            
            %get out other signals
            current_H2B = H2B_store(track,:);
            FRET = YFP_store(track,:)./CFP_store(track,:);
            H2B_norm = current_H2B/max(current_H2B);
            
            %filter out the poorly expressing cells
            if max(FRET) <= resting_FRET_condition+1 && min(FRET) >= resting_FRET_condition-1
                
                %ensures that the mitosis event occurs, and it happens
                %after 10 frames into the movie and before 5 frames from
                %the end
                if ~isempty(mitoses) && mitoses(1) > 10 && mitoses(end) < 170
                    
                    %if mitosis events happen multiple times, we are
                    %cycling through them
                    for k = 1:length(mitoses)
                        
                        %gets out the mean "resting FRET", or the average FRET Ratio 5 frames before
                        %mitosis (~2h to 1h)
                        resting_FRET = [resting_FRET mean(FRET(mitoses(k)-9:mitoses(k)-5))];
                        
                        %finds the maximum value between mitosis and 1h
                        %after mitosis
                        peak_FRET = [peak_FRET max(FRET(mitoses(k):mitoses(k)+5))];
                        
                    end
                end
                
            end
        end
        
    end
    
    %calculates zinc spike from the peak and resting FRET
    zn_spike = peak_FRET-resting_FRET;
    
    %uses write_spike_spreadsheet function to write the spike information
    %to a spreadsheet, using the column number as the media type
    
    condition = condition_cell{c};
    write_spike_spreadsheet(resting_FRET,peak_FRET,zn_spike,condition,sensor,output_str)
    
end