%Contents of this code
%As of 09/22/2021

%main_file
%This is where you put your experiment-specific information such as the
%date of the experiment, sensor type, cell-specific information such as
%KDs, colors of plotting each treatment

%step0_cellproliferation plots the average cell count at each frame over
%time and plots each condition in a different color

%step0_fig1C_MTF1KD plots cells and the mitosis events of each cell as red
%dots in addition to the FRET ratio of the cells over time

%step1_resting_FRET_all_treatments_asynchronous averages the FRET ratios of
%all tracks as they asynchronously cycle. The law of large numbers states
%this will average the tracks to the "theoretical mean" of the FRET ratio
%in each  condition. This serves as a downstream proxy for "resting FRET"
%that the pipeline uses for downstream calculations and filters.

%step2_mean_resting_FRET_all_treatments_align takes all FRET ratios and
%aligns them to the first  mitosis the tracks experience. It then averages
%these tracks and plots them on an axis labelled "hours relative to
%mitosis"

%step3_all_tracks_all_treatments_aligned_first_mit does the same as above,
%but does not average and plots each individual track in an m x n subplot