%This function writes zinc spike data to a spreadsheet, where the condition
%in variable 'condition' corresponds to a separate sheet in the spreadsheet

function write_spike_spreadsheet(resting_FRET,peak_FRET,zn_spike,condition,sensor,outdir)

resting_FRET = resting_FRET';
peak_FRET = peak_FRET';
zn_spike = zn_spike';

all_zn_data = [resting_FRET peak_FRET zn_spike];

writematrix(all_zn_data, [outdir,'/', 'zn_spike_', sensor, '.xls'], 'Sheet', condition)

end