function proliferation_vec = get_proliferation(all_signals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%initializing the storage vector
existCell_well = [];

if ~isempty(all_signals)
    
    well = all_signals{1};
    %number of tracks in each well
    tracks = length(well);
    
    existVec = [];
    
    for i=1:tracks
        %cycle through tracks
        
        track = well{i}.ellipse_id';
        %index the ellipse id, which is a proxy for when the ellipse
        %exists, therefore when the cell exists
        existVec=[existVec;track];
        %store this value in existVec vertically
        
    end
    
    existCell_well = [existCell_well;existVec];
    %store this in the cell for the entire well vertically
end



%put each well into a cell, where the row corresponds to the
%treatment

cell_exists = ~isnan(existCell_well);
%get a logical T/F where the vector is or isn't a NaN

cell_count = sum(cell_exists, 1);

proliferation_vec = cell_count / cell_count(1);

end

