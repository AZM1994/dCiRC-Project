%%% Load 'CIRC_Data_Cell.mat'
%%% Select certain data entry for specific genotype
%%% Detrend, smooth, and normalize the data entry
%%% Prepare data for the optimization algorithms
%%% Written by Zheming An
clear, close all

%% Load race tube data from 'CIRC_Data_Cell.mat'
CIRC_Data_Cell = load('CIRC_Data_Cell.mat');
CIRC_Data_Cell = CIRC_Data_Cell.CIRC_Data_Cell;
Size_CIRC_Data_Cell = size(CIRC_Data_Cell,1);
data_range = 100; % check 'CIRC_Data_Cell.mat', select the data entries for specific genotypes,
                  % e.g., data_range = 100 corresponds to genotype JW18_2,
                  % data_range = 99:101 corresponds to the average of three
                  % entries of genotype JW18

allCIRC_Data = 0; 

for i = data_range
CIRC_LD_Phase_Data_1  = CIRC_Data_Cell{i,3}(1:end);
CIRC_LD_Time_Data     = CIRC_Data_Cell{i,2}(1:end);

figure(1);
plot(CIRC_LD_Time_Data, CIRC_LD_Phase_Data_1, 'LineWidth',2)
xlabel('time/hours')
ylabel('rhythmic data')
hold on

title_string = sprintf('Data in row(s) %d - %d (%s)', data_range(1), data_range(end), char(CIRC_Data_Cell{data_range(1),1})); 
title(title_string)
set(gca, 'FontSize', 24) 
allCIRC_Data = allCIRC_Data + CIRC_LD_Phase_Data_1; 
end

% Data truncation (based on the features of each data entry, remove the front and end of the data)
data_range_truncation = 15:685; % Choose manually to make the data processable
allCIRC_Data = allCIRC_Data(data_range_truncation);
CIRC_LD_Time_Data = CIRC_LD_Time_Data(data_range_truncation);

% Taking average if there are more than one entries selected
CIRC_LD_Phase_Data = allCIRC_Data / length(data_range);
figure(1);
plot(CIRC_LD_Time_Data, CIRC_LD_Phase_Data, 'k', 'LineWidth',3)



%% Detrend the selected race tube data by subtracting linear/quadratic/cubic trend
Detrend_LD_Data = detrend(CIRC_LD_Phase_Data,3); % Detrand by subtracting the cubic trend
figure(2);
plot(CIRC_LD_Time_Data, CIRC_LD_Phase_Data,'b',...
     CIRC_LD_Time_Data, CIRC_LD_Phase_Data - Detrend_LD_Data, '--k',...
     CIRC_LD_Time_Data, Detrend_LD_Data, 'r',...
     'LineWidth',3)
xlabel('time/hours')
ylabel('rhythmic data')
legend('Original Data', 'Trend', 'Detrended Data')
set(gca, 'FontSize', 24)



%% Smooth the detrended race tube data with function csaps (Cubic smoothing spline)
figure(3);
plot(CIRC_LD_Time_Data, Detrend_LD_Data,'r', 'LineWidth',3)
hold on

factor = 0.005; % factor is called the smoothing parameter
s0 = csaps(CIRC_LD_Time_Data, Detrend_LD_Data,factor);
fnplt(s0,'b',3); % plot the smoothed data
Smoothed_LD = fnval(s0,CIRC_LD_Time_Data);
xlabel('time/hours')
ylabel('rhythmic data')
legend('Detrended Data', 'Smoothed Data')
set(gca, 'FontSize', 24)



%% Normalize the smoothed race tube data to make the amplitude equals to 1
fh = figure(4); 
fh.WindowState = 'maximized';
plot(CIRC_LD_Time_Data, Detrend_LD_Data,'g', 'LineWidth',3)
hold on
xticks(0:24:120)
xlabel('time/hours')
ylabel('rhythmic data')
set(gca,'FontSize',24)
ylim([-1.1 1.1])
title_string_2 = sprintf('Data in %d - %d (%s), range in %d - %d,\n smoothing parameter = %d', ...
    data_range(1), data_range(end), char(CIRC_Data_Cell{data_range(1),1}),...
    data_range_truncation(1), data_range_truncation(end), factor); 
title(title_string_2)

% Find peaks of smoothed data
MinPeakDistance_pos = 100;
MinPeakDistance_neg = 100;
poly_fit_degree_pos = 5;
poly_fit_degree_neg = 5;

[pos_peak_value, pos_peak_time_index] = findpeaks(Smoothed_LD,'MinPeakDistance',MinPeakDistance_pos); 
[neg_peak_value, neg_peak_time_index] = findpeaks(-Smoothed_LD,'MinPeakDistance',MinPeakDistance_neg); 
pos_peak_time = CIRC_LD_Time_Data(pos_peak_time_index); 
neg_peak_time = CIRC_LD_Time_Data(neg_peak_time_index); 
number_of_peaks = min(length(pos_peak_time), length(neg_peak_time)); 
number_of_peaks_pos = length(pos_peak_time); 
number_of_peaks_neg = length(neg_peak_time); 

plot(CIRC_LD_Time_Data, Smoothed_LD,'--k', 'LineWidth',3)        

% generate the polynomial fitting line of the positive and negative peaks
amplitude_fitting_pos = polyfit([CIRC_LD_Time_Data(1);pos_peak_time(1 : number_of_peaks_pos);CIRC_LD_Time_Data(end)],...
    [pos_peak_value(1);pos_peak_value(1 : number_of_peaks_pos);pos_peak_value(end)],poly_fit_degree_pos); 
amplitude_fitting_neg = polyfit([CIRC_LD_Time_Data(1);neg_peak_time(1 : number_of_peaks_neg);CIRC_LD_Time_Data(end)],...
    -[neg_peak_value(1);neg_peak_value(1 : number_of_peaks_neg);neg_peak_value(end)],poly_fit_degree_neg); 
amplitude_regression_pos = polyval(amplitude_fitting_pos, CIRC_LD_Time_Data); 
amplitude_regression_neg = polyval(amplitude_fitting_neg, CIRC_LD_Time_Data); 
plot(CIRC_LD_Time_Data, amplitude_regression_pos, '-.r'); 
plot(CIRC_LD_Time_Data, amplitude_regression_neg, '-.b'); 

% Normalize the data with polynomial fitting lines
Normalized_amplitude_LD = zeros(length(Smoothed_LD),1); 

for i = 1 : length(Smoothed_LD)
    if Smoothed_LD(i) > 0
        Normalized_amplitude_LD(i) = Smoothed_LD(i) ./ amplitude_regression_pos(i);
    else
        Normalized_amplitude_LD(i) = Smoothed_LD(i) ./ -amplitude_regression_neg(i); 
    end
end

% Find peaks of the normalized data
[pos_peak_value_amplitude, pos_peak_time_index_amplitude]= findpeaks(Normalized_amplitude_LD,'MinPeakDistance',MinPeakDistance_pos); 
[neg_peak_value_amplitude, neg_peak_time_index_amplitude]= findpeaks(-Normalized_amplitude_LD ,'MinPeakDistance',MinPeakDistance_neg); 
pos_peak_time_amplitude = CIRC_LD_Time_Data(pos_peak_time_index_amplitude); 
neg_peak_time_amplitude = CIRC_LD_Time_Data(neg_peak_time_index_amplitude); 
number_of_peaks_amplitude = min(length(pos_peak_time_amplitude), length(neg_peak_time_amplitude)); 
Normalized_amplitude_LD(Normalized_amplitude_LD > ...
    max(pos_peak_value_amplitude)) = NaN; 
Normalized_amplitude_LD(Normalized_amplitude_LD < ...
    -max(neg_peak_value_amplitude)) = NaN; 

max_NA = max(Normalized_amplitude_LD); 
min_NA = min(Normalized_amplitude_LD); 
Normalized_amplitude_LD = Normalized_amplitude_LD * (max_NA - min_NA) / 2 + (max_NA + min_NA) / 2; 
Normalized_amplitude_LD = (Normalized_amplitude_LD - (max_NA + min_NA) * (max_NA - min_NA) / 4 - (max_NA + min_NA) / 2) / ((max_NA - min_NA)^2 / 4); 
plot(CIRC_LD_Time_Data, Normalized_amplitude_LD, 'm', 'LineWidth', 5)

legend('Detrended Data', 'Smoothed Data', 'Amplitude Regression (pos)', 'Amplitude Regression (neg)',...
       'Normalized Data', 'Location', 'southeast', 'AutoUpdate','off')
plot(CIRC_LD_Time_Data, zeros(length(CIRC_LD_Time_Data),1), '--', 'Color',[0.5 0.5 0.5], 'LineWidth',3)
plot(CIRC_LD_Time_Data, zeros(length(CIRC_LD_Time_Data),1) + 1,'--k', 'LineWidth',2)
plot(CIRC_LD_Time_Data, zeros(length(CIRC_LD_Time_Data),1) - 1,'--k', 'LineWidth',2)
scatter(pos_peak_time(1 : number_of_peaks_pos), pos_peak_value(1 : number_of_peaks_pos), 200, 'filled', 'r');
scatter(neg_peak_time(1 : number_of_peaks_neg), -neg_peak_value(1 : number_of_peaks_neg), 200, 'filled', 'b'); 
scatter(pos_peak_time_amplitude(1 : number_of_peaks_amplitude), pos_peak_value_amplitude(1 : number_of_peaks_amplitude), 100, 'filled', 'r');
scatter(neg_peak_time_amplitude(1 : number_of_peaks_amplitude), -neg_peak_value_amplitude(1 : number_of_peaks_amplitude), 100, 'filled', 'b'); 

% save normalized race tube data for optimization process
normalized_data_saved = [CIRC_LD_Time_Data(~isnan(Normalized_amplitude_LD)),...
    Normalized_amplitude_LD(~isnan(Normalized_amplitude_LD))]; 
save('normalized_data_saved','normalized_data_saved')
