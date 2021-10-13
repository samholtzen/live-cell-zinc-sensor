%% 

clearvars
[filename,path] = uigetfile();
S = load([path,filename]);
fieldname = fieldnames(S);
signals = S.(fieldname{1});
%% 

%Date of experiment
experiment_date = '20170120';

%Sensor used
sensor = 'nls';

%Does it bind zinc? Or is it dead? (You can also put any other info here,
%like knockdowns or other cell-type specific identifiers)
binding = 'dead';

%What are the conditions in each of your columns (list in order)?
condition_cell = {'ZD','ZD + 2µM Zn','Minimal Media','ZR30','ZR50', 'Non-chelexed horse serum'};
colors_cell = {'#B47EB3','#083D77','#FFD5FF','#92D1C3','#8BB8A8','#5D4E6D'};

%How many columns are you plotting
column_number = length(condition_cell);

%How many replicates have you done in each condition?
replicate_number = size(signals,1);

%which conditions do you want to plot?
conditions_to_plot = [1, 2, 3, 4, 5, 6];

%is there a media change? if so, what frame did it occur?

media_frame = 1;

%% 

output_str = [sensor,'_',binding,'_','output'];
if ~exist(output_str,'dir')
    mkdir(output_str)
end



%Run step0_cellproliferation

step0_cellproliferation
saveas(gcf,[output_str,'/','cell_proliferation', experiment_date, sensor, binding, '.png'])
close(gcf)

disp('Finished step0_cellproliferation and wrote file to output folder.')

%Run step0_fig1C_MTFKD
 
step0_fig1C_MTF1KD
saveas(gcf,[output_str,'/','mitosis_red_dots', experiment_date, sensor, binding, '.png'])
close(gcf)

disp('Finished step0_fig1C_MTF1KD and wrote file to output folder.')

%Run step1_resting_FRET_all_treatments_asynchronous
step1_resting_FRET_all_treatments_asynchronous
saveas(gcf,[output_str,'/','resting_fret_asynchronous', experiment_date, sensor, binding, '.png'])
close(gcf)

disp('Finished step1_resting_FRET_all_treatments_asynchronous and wrote file to output folder.')

%Run step2_resting_FRET_all_treatments_align
step2_mean_resting_FRET_all_treatments_align
saveas(gcf,[output_str,'/','resting_fret_aligned_mitosis', experiment_date, sensor, binding, '.png'])
close(gcf)

disp('Finished step2_mean_resting_FRET_all_treatments_align and wrote file to output folder.')

%Run step3_all_tracks_all_treatments_aligned_first_mit
step3_all_tracks_all_treatments_aligned_first_mit
saveas(gcf,[output_str,'/','all_tracks_fret_aligned_mitosis', experiment_date, sensor, binding '.png'])
close(gcf)

disp('Finished step3_all_tracks_all_treatments_aligned_first_mit and wrote file to output folder.')

%Run step4_zinc_spike_spreadsheet
step4_zinc_spike_spreadsheet

disp('Finished step4_zinc_spike_spreadsheet and wrote file to output folder.')