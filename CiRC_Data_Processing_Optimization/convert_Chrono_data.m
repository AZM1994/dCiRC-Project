%%% convert Chrono race tube data 'P4414.txt' to 'CIRC_Data_Cell.mat'
%%% 'CIRC_Data_Cell.mat' is a 411 by 3 cell
%%% First column contains the name of the genotypes, 1-3 entries for each.
%%% Second column contains the time data of each experiment
%%% Third column contains the phase data of each experiment
%%% Written by Zheming An
clear, close all

%% Generate CIRC_Data_Cell.mat
file_name    = 'P4414.txt';
CIRC_LD_Data = readmatrix(file_name);
opts         = detectImportOptions(file_name,'Range','1:3');
Name_Row     = opts.VariableNames;

CIRC_Data_Cell      = cell(floor(size(CIRC_LD_Data,2) / 4), 3);

for i = 1 : size(CIRC_Data_Cell,1)
data_column_index   = i * 4 - 2;
CIRC_LD_Phase_Data  = (CIRC_LD_Data(1:end,data_column_index));
CIRC_LD_Time_Data   = (CIRC_LD_Data(1:end,data_column_index - 1) / 1440 * 24);
CIRC_Data_Cell{i,1} = Name_Row(1, i * 4 - 3);
CIRC_Data_Cell{i,2} = CIRC_LD_Time_Data;
CIRC_Data_Cell{i,3} = CIRC_LD_Phase_Data;
end

save('CIRC_Data_Cell.mat', 'CIRC_Data_Cell');