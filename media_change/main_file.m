%% Selecting files to analyze
clearvars
[filename,path] = uigetfile('MultiSelect','on');

if ~iscell(filename)
    % Workaround for the load function 
    filename = {filename};
end



for i=1:length(filename)
    
    S = load([path,filename{i}]);
    fieldname = fieldnames(S);
    struct_cell = S.(fieldname{1});
    
    rng(1)
    split_filename = split(filename{i},'_');
    %Date of experiment
    experiment_date = split_filename{1};
    
    %Sensor used
    sensor = 'nes';
    
    %Does it bind zinc? Or is it dead? (You can also put any other info here,
    %like knockdowns or other cell-type specific identifiers)
    colors_cell = {'#33637a','#3DA5D9','#0C0808','#CE9600','#ED8A07','#B35C05'};
    cell_type = split_filename{2};
    
    switch cell_type
        
        case {'scr','all'}
            condition_cell = {'ZD3','ZD2','Minimal Media','ZR15','ZR30','ZR50'};
            conditions_to_plot = [1,2,3,4,5,6];
            
            
        case 'kd'
            condition_cell = {'ZD3','ZD2','Minimal Media','ZR15','ZR30','ZR50'};
            conditions_to_plot = [1, 2, 3, 4, 5,6];
            
    end
    
    switch experiment_date
        
        case '20220713'
            media_change = 118;
            
        case '20220722'
            media_change = 127;
            
        otherwise
            media_change = 120;
            
    end
    
    output_str = [sensor,'_',cell_type,'_','output'];

    
    
    % Run parameters for filtering YFP signal. These are empirically tested to
    % ensure we remove cells that aren't expressing the construct.
    
    FRET_min = 3;
    FRET_max = 8;
    
    y_min = 3.5;
    y_max = 7;
    
%% Uncomment below to export graphs  

%     if ~exist(output_str,'dir')
%         mkdir(output_str)
%     end

%     step0_mit_red_dots
%     savefig(gcf,[output_str,'/',experiment_date,'_mitosis_red_dots_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_mitosis_red_dots_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     %Run step1_resting_FRET_all_treatments_asynchronous
%     step1_resting_FRET_all_treatments_asynchronous
%     savefig(gcf,[output_str,'/',experiment_date,'_resting_fret_asynchronous_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_resting_fret_asynchronous_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     step2e_zinc_estimation
%     savefig(gcf,[output_str,'/',experiment_date,'_estimated_zinc_asynchronous_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_estimated_zinc_asynchronous_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     %Run step2_resting_FRET_all_treatments_align
%     step2_mean_resting_FRET_all_treatments_align
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     savefig(gcf,[output_str,'/',experiment_date,'_resting_fret_align_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_resting_fret_align_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     %Run step4_zinc_spike_spreadsheet
%     step4_media_change_spreadsheet
%     
%     step5_additional_media_change
    
end