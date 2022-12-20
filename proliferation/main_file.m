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
    split_filename = split(filename{i},'_');
    
    %Date of experiment
    experiment_date = split_filename{1};
    
    %Sensor used
    sensor = 'nes';
    
    cell_type = split_filename{2};
    
    colors_cell = {'#33637a','#3DA5D9','#0C0808','#CE9600','#ED8A07','#B35C05'};
    
    switch cell_type
        
        case {'scr','all'}
            condition_cell = {'ZD3','ZD2','Minimal Media','ZR15','ZR30','ZR50'};
            conditions_to_plot = [1,2,3,4,5,6];
            y_min = 4;
            y_max = 6.5;
            ifcrop = true;
            
        case 'kd'
            condition_cell = {'ZD3','ZD2','Minimal Media','ZR15','ZR30','ZR50'};
            conditions_to_plot = [1, 2, 3, 4, 5,6];
            y_min = 4;
            y_max = 6.5;
            ifcrop = true;
            
        case {'dead','wt'}
            colors_cell = {'#3DA5D9','#3DA5D9','#0C0808','#CE9600','#ED8A07','#B35C05'};
            condition_cell = {'ZD2','ZD + 3µM Zn','Minimal Media','ZR15','ZR30','No Chx'};
            conditions_to_plot = [1,3, 4, 5];
            y_min = 3;
            y_max = 6.5;
            ifcrop = false;
            
    end
    
    
    output_str = [sensor,'_',cell_type,'_','output'];
    
    
    
    % Run parameters for filtering YFP signal. These are empirically tested to
    % ensure we remove cells that aren't expressing the construct.
    
    FRET_min = 3;
    FRET_max = 8;
    
    %% Uncomment to export figures
    
%     if ~exist(output_str,'dir')
%         mkdir(output_str)
%     end
%     
%     
%     step0_cellproliferation
%     savefig(gcf,[output_str,'/',experiment_date,'_proliferation_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_proliferation_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     step0_fig1C_MTF1KD
%     savefig(gcf,[output_str,'/',experiment_date,'_mitosis_red_dots_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_mitosis_red_dots_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     step1_resting_FRET_all_treatments_asynchronous
%     savefig(gcf,[output_str,'/',experiment_date,'_resting_fret_asynchronous_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_resting_fret_asynchronous_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     step2_mean_resting_FRET_all_treatments_align
%     savefig(gcf,[output_str,'/',experiment_date,'_resting_fret_aligned_mitosis_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_resting_fret_aligned_mitosis_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     step2a_mean_aligned_skinny
%     savefig(gcf,[output_str,'/',experiment_date,'_resting_fret_aligned_skinny_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_resting_fret_aligned_skinny_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     step2b_subsample_mean_align
%     savefig(gcf,[output_str,'/',experiment_date,'_subsample_mean_align_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_subsample_mean_align_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     step2c_mean_align_resting_dash
%     savefig(gcf,[output_str,'/',experiment_date,'_align_resting_dash_', sensor,'_', cell_type, '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_align_resting_dash_', sensor,'_', cell_type, '.png'],'Resolution',300)
%     close(gcf)
%     
%     step3_all_tracks_all_treatments_aligned_first_mit
%     savefig(gcf,[output_str,'/',experiment_date,'_all_tracks_fret_aligned_mitosis_', sensor,'_', cell_type '.fig'])
%     exportgraphics(gcf,[output_str,'/',experiment_date,'_all_tracks_fret_aligned_mitosis_', sensor,'_', cell_type '.png'],'Resolution',300)
%     close(gcf)
%     
%     step4_zinc_spike_spreadsheet
    
end