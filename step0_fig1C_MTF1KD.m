%This plots mitosis events as red dots over traces of FRET ratio in dark
%gray

%columns are different media conditions, rows are replicates

for c=1:column_number
    %loop through the columns (conditions)
    
    
    for nd = 1:replicate_number
        %loop through the replicates and extract the FRET values, H2B
        %values, and mitosis values
        
        wells = signals(nd,c,1);
        
        CFP_store = [];
        YFP_store = [];
        H2B_store = [];
        mitosis_store = [];
        
        %% Find and join mother/daughter pairs to create full lineage
        
        for j=1:numel(wells)
            
            %cycling through wells
            
            well = wells{j};
            
            [all_full_traces_CFP, all_full_traces_YFP, all_full_traces_H2B, all_mitosis,num_frames] = link_mother_daughter(well,sensor);
            
            CFP_store = [CFP_store;all_full_traces_CFP];
            YFP_store = [YFP_store;all_full_traces_YFP];
            H2B_store = [H2B_store;all_full_traces_CFP];
            mitosis_store = [mitosis_store; all_mitosis];
        end
        %% Get mitosis events and use those to plot FRET ratio and red dots
        
        %create random tracks to plot
        random = randi([1,length(mitosis_store)],1,30);
        
        %initialize storage vectors
        first_mitoses = [];
        first_mit_FRET = [];
        
        %loop through the cells in one condition
        for i = 1:length(mitosis_store)
            
            %find mitosis events
            mitoses = mitosis_store{i};
            
            %convert frames to hours
            frame_vec = 1:num_frames;
            
            %calculate FRET Ratio (YFP/CFP)
            current_FRET = YFP_store(i,:)./CFP_store(i,:);
            
            if max(current_FRET) <= 8 && min(current_FRET) >= 3
                %filter out aberrant sensor expression by constraining to a
                %reasonable FRET value in the DR of the sensor
                
                
                if  ~isempty(mitoses) && ismember(i,random)
                    %if elliptrack recorded a mitosis event
                    subplot(2,3,c)
                    h(1) = plot(frame_vec./5, current_FRET, 'Color', [0.25 0.25 0.25]);
                    hold on
                    
                    % store the mitosis events and associated FRET in a vector to plot at the
                    % end
                    first_mitoses = [first_mitoses; mitoses(1)];
                    first_mit_FRET = [first_mit_FRET; current_FRET(mitoses(1))];
                    
                elseif ismember(i,random)
                    %elliptrack didn't see a mitosis event, plot anyway
                    subplot(2,3,c)
                    h(1) = plot(frame_vec./5, current_FRET, 'Color', [0.25 0.25 0.25]);
                    
                    hold on
                    
                end
                
            end
            
            
        end
    end
    
    %plot the stored mitosis events
    subplot(2,3,c)
    h(2) = plot(first_mitoses./5,first_mit_FRET, '.r','MarkerSize',10);
    hold on
    
end


for condition = 1:column_number
        subplot(2,3,condition)
        xlabel('Time (hours)')
        ylabel('FRET Ratio')
        title(['Mitosis Events,',condition_cell{condition}])
end