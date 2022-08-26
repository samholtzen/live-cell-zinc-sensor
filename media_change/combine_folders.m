% Write a script that loops through all the output files and puts each well
% into a cell array like normal

num_rows = 5;
num_cols = 12;

all_signals_temp = cell(num_rows, num_cols);

for i = 1:num_rows
    
    for j= 1:num_cols
        
        filestem = ['results/', num2str(i), '_', num2str(j),'_1/'];
        
        if isequal(exist([filestem,'signals.mat'], 'file'),2)
            
            load([filestem,'signals.mat'])
            
            all_signals_temp{i,j} = all_signals{i,j};
            
            
            
        end
        
        disp(['Added element in row',i,', column',j])
        
    end
    
end

all_signals = all_signals_temp;

%save('signals.mat','all_signals','-v7.3')