%figure 1B on printout of figures for MTF1 KD paper

%initializing the storage vectors and cells
existCell = {};
existCell_well = [];
existVecLog_temp = [];
existVecLog = [];
frames = 175;

for c=1:column_number
    
    %initializing the storage vector
    existCell = {};
    existCell_well = [];
    %looping through eac well of the column
    for nd = 1:replicate_number
        
        %pulling cell information from each cell of the storage array,
        %where nd is the row number, and c is the treatment condition
        wells = signals{nd,c,1};
        
        %number of tracks in each well
        tracks = length(wells);
        
        existVec = [];
        
        for i=1:tracks
            %cycle through tracks
            
            track = wells{i}.ellipse_id';
            %index the ellipse id, which is a proxy for when the ellipse
            %exists, therefore when the cell exists
            existVec=[existVec;track];
            %store this value in existVec vertically
            
        end
        
        existCell_well = [existCell_well;existVec];
        %store this in the cell for the entire well vertically
        existCell = [existCell;{existCell_well}];
        
        
    end
    
    %put each well into a cell, where the row corresponds to the
    %treatment
    existVecLog = [existVecLog,existCell];
end

existVecLog = cellfun(@isnan,existVecLog,'UniformOutput',false);
%get a logical T/F where the vector is or isn't a NaN
frameVec = (1:frames)/5;

for condition = 1:column_number
    for rep = 1:replicate_number
        cell_count = smooth(sum(~existVecLog{rep,condition}));
        norm_cell_count = cell_count/cell_count(1);
        existVecLog(rep,condition) = {norm_cell_count};
        %Getting the info out of the "existVecLog" variable and smoothing so the
        %graph looks okay
    end
    
    cell_count_plot = [existVecLog{:,condition}];
    cell_count_plot_mean = mean(cell_count_plot,2);
    plot(frameVec,cell_count_plot_mean,'Color',colors_cell{condition},'LineWidth',2)
    hold on
end

legend(condition_cell,'Location','Northwest')
ylabel('Normalized cell count')
xlabel('Time (Hours)')


