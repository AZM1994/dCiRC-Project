%% For genotype JW59, Generate phase response curve (PRC) with the best fitting dCiRC paramters

clear;close all;

global plot_time T_final light_profiles_LD light_profiles_DD light_profiles_index zeit_t organism_names ...
       light_period light_off light_on save_analytical_phase_difference_LD save_analytical_phase_difference_DD
%% Load data
% best fitting dCiRC parameters for all genotypes
optimized_x_matrix = [21.0	0.6477	1.5737	1.40	0.05	-0.27	18.00	96
                    21.0	0.2436	0.7246	0.00	0.05	-0.24	18.00	96
                    21.0	0.3696	1.5247	0.05	0.10	-0.0321	18.00	96
                    21.2	0.3548	1.9021	0.00	0.10	-0.22	18.00	96
                    21.3	1.4254	1.7593	1.05	0.05	-0.17	18.00	96
                    21.3	1.0021	1.9621	0.00	0.05	-0.15	18.00	96
                    21.5	0.2958	0.6504	0.05	0.05	0.02	18.00	96
                    21.7	0.2658	0.4019	0.20	0.10	0.06	18.00	96
                    21.7	0.3373	0.3216	1.00	0.05	-0.28	18.00	96
                    21.7	0.1559	0.1701	0.85	0.60	-0.03	18.00	96
                    21.8	0.1536	0.1060	1.40	0.40	-0.06	18.00	96
                    22.0	0.5679	0.5591	0.30	0.10	-0.15	18.00	96
                    22.2	0.9090	1.0826	0.40	0.25	-0.07	18.00	96
                    22.2	0.3997	1.1089	0.10	0.25	0.08	18.00	96
                    22.2	0.3784	1.9092	0.00	0.57	0.01	18.00	96
                    22.2	0.4937	1.9138	0.05	0.15	0.41	18.00	96
                    22.3	0.3879	1.0788	0.20	0.40	-0.20	18.00	96
                    22.5	0.3429	0.7786	0.15	0.20	0.00	18.00	96
                    22.7	0.6499	1.9883	0.00	0.10	0.33	18.00	96
                    22.7	0.3258	0.7034	0.26	0.35	0.27	18.00	96
                    23.2	0.2915	0.9880	0.20	0.10	-0.02	18.00	96
                    23.5	0.5339	0.6208	0.60	0.30	0.17	18.00	96
                    23.8	0.7301	1.2864	0.00	0.05	0.24	18.00	96
                    23.8	0.4537	0.8656	0.60	0.20	0.28	18.00	96
                    22.0	0.3636	1.1602	1.00	0.14	0.00	18.00	96
                    30.0	0.0394	0.0056	1.40	1.85	0.23	18.00	96];

optimized_x = optimized_x_matrix;
organism_names  = ["JW168", "JW224", "JW60",  "JW260", "D117",...
                  "D119",  "D116",  "JW18",  "JW172", "P4483",...
                  "JW220", "JW200", "JW162", "JW176", "JW238",...
                  "JW261", "JW180", "JW169", "JW228", "P4463",...
                  "JW24",  "JW59",  "JW161", "P4469", "JW22",...
                  "DBP338"];
organism_index  = 22;
Genotype        = organism_names(organism_index);

%% parameters
T_initial       = 0;
T_interval      = 0.001;
T_final         = optimized_x(organism_index,8);
light_period    = 840; 
light_off       = 24; 
light_on        = 0;
number_of_timesteps = (T_final - T_initial) / T_interval;
plot_time       = T_initial : T_interval : (T_final - T_interval);
zeit_t          = 24;

%% mimic PRC light pulse settings
seq_01 = 1:24;
seq_02 = repmat(0.5,[1,length(seq_01)]);
light_profiles_hourly_LD = [seq_02',seq_01'];
light_profiles_LD  = light_profiles_hourly_LD * 2 * pi / zeit_t; % dawn dusk

seq_03 = zeros([1,length(seq_01)]);
light_profiles_hourly_DD = [seq_03',seq_01'];
light_profiles_DD  = light_profiles_hourly_DD * 2 * pi / zeit_t; % dawn dusk

save_analytical_phase_difference_LD = zeros(T_final/zeit_t,size(light_profiles_LD,1));
save_analytical_phase_difference_DD = zeros(T_final/zeit_t,size(light_profiles_DD,1));

%% simulate the phase changes 
for light_profiles_index = 1 : size(light_profiles_LD,1)
    CIRC_LD_Cost_for_plot_best_fitting_LD(optimized_x,organism_index,T_interval,number_of_timesteps)
end

for light_profiles_index = 1 : size(light_profiles_DD,1)
    CIRC_LD_Cost_for_plot_best_fitting_DD(optimized_x,organism_index,T_interval,number_of_timesteps)
end

%% Plot PRC
fh = figure(1);
fh.WindowState = 'maximized';

% shift the simulated phase to generate PRC (align the time axis)
record_after_pulse = 2;
phase_shift = save_analytical_phase_difference_DD(record_after_pulse,:) - save_analytical_phase_difference_LD(record_after_pulse,:);

for i = 1 : length(phase_shift)
    if light_profiles_hourly_DD(i,2) < 6
        light_profiles_hourly_DD(i,2) = light_profiles_hourly_DD(i,2) + zeit_t;
    end
end

PRC_data = [light_profiles_hourly_DD(:,2),phase_shift'];
PRC_data_new = sortrows(PRC_data);
PRC_time = PRC_data_new(:,1);
PRC_phase = PRC_data_new(:,2);
PRC_time = PRC_time * zeit_t / optimized_x_matrix(1,1);
plot([PRC_time;PRC_time(end)+zeit_t / optimized_x_matrix(1,1)], [PRC_phase;PRC_phase(1)],'k','LineWidth',10)
hold on
plot(PRC_time,zeros(1, length(PRC_time)),'--k','LineWidth',3)
text(29,0.35,'advance','Color','r','FontSize',42,'Rotation',0);
text(20.5,-1.4,'delay','Color','b','FontSize',42,'Rotation',0);
xticks(PRC_time(1):3*zeit_t / optimized_x_matrix(1,1):PRC_time(end)+zeit_t / optimized_x_matrix(1,1))
xticklabels({'CT0','CT3','CT6','CT9','CT12','CT15','CT18','CT21','CT0'})
xlim([PRC_time(1) PRC_time(end)+zeit_t / optimized_x_matrix(1,1)])
ylim([-4 2])
xlabel('circadian time in $\tau$/24 hour','Interpreter','latex')
ylabel('phase shift/hour','Interpreter','latex')
legend(['PRC (' + Genotype + ')'],'Location','best')
set(gca, 'FontSize',42)
set(gca,'LineWidth',2)



%% dCiRC model
function CIRC_LD_Cost_for_plot_best_fitting_LD(optimized_x,organism_index,T_interval,number_of_timesteps)
global T_final light_profiles_LD light_profiles_index zeit_t ...
       light_period light_off light_on save_analytical_phase_difference_LD

% preallocate
ALL_clocks                  = zeros(number_of_timesteps,2);
zeitgeber_clock             = zeros(number_of_timesteps,1);
Normalized_theta_Z          = zeros(number_of_timesteps,1);
Normalized_theta_I          = zeros(number_of_timesteps,1);
CIRC                        = zeros(number_of_timesteps,1);
light                       = zeros(number_of_timesteps,1);
omega_I                     = zeros(number_of_timesteps,1);
domega_I                    = zeros(number_of_timesteps,1);

% Parameters
tau_DD          = optimized_x(organism_index,1);
C_ZG            = optimized_x(organism_index,2);
elastic_factor  = optimized_x(organism_index,3);
s               = optimized_x(organism_index,4);
a               = optimized_x(organism_index,5);
phi_0           = optimized_x(organism_index,6);
zeit_shift      = optimized_x(organism_index,7);
omega_DD        = 2 * pi ./ tau_DD;
initial_theta_I = phi_0;
initial_omega_I = omega_DD;
CIRC_shift      = 0;

% light profiles
epsilon         = light_profiles_LD(light_profiles_index,1); % radius of the entrainment window
light_middle    = light_profiles_LD(light_profiles_index,2); % center of the entrainment window

%% Simulations
theta_I_initial = initial_theta_I; 
omega_I_initial = initial_omega_I; 
initial_conditions = [theta_I_initial, omega_I_initial];
ALL_clocks(1,:) = initial_conditions;

for i = 1:number_of_timesteps - 1
    zeitgeber_clock(i) = (i-1-(18-zeit_shift)/T_interval)*2*pi/(zeit_t/T_interval);
    k1= RK4_circadian_CIRC(ALL_clocks(i,:));
    k2_clocks = ALL_clocks(i,:) + 0.5*T_interval*k1;

    k2= RK4_circadian_CIRC(k2_clocks);
    k3_clocks = ALL_clocks(i,:) + 0.5*T_interval*k2;

    k3= RK4_circadian_CIRC(k3_clocks);
    k4_clocks = ALL_clocks(i,:) + T_interval*k3;

    k4= RK4_circadian_CIRC(k4_clocks);
    ALL_clocks(i+1,:) =  ALL_clocks(i,:) + T_interval*(k1 + 2*k2 + 2*k3 + k4) / 6.0;             
end
    analytical_zeitgeber_signal = sin(zeitgeber_clock);
    theta_I_signal = sin(ALL_clocks(:,1));

% find the time values corresponding to the peaks of the signal
    [~, n_time]= findpeaks(theta_I_signal,'MinPeakHeight',0.99);
    [~, analytical_z_time]= findpeaks(analytical_zeitgeber_signal,'MinPeakHeight',0.99);

    analytical_days = size(analytical_z_time,1);
    analytical_phase_difference = zeros(analytical_days,1);
    
    for j = 1:analytical_days
        [~, analytical_phase_index] = min(abs(analytical_z_time(j) - n_time));
        analytical_phase_difference(j) = (n_time(analytical_phase_index) - analytical_z_time(j)) .* T_interval;
    end
    
    entrained_time = T_final;
    average_POE = mean(analytical_phase_difference);
    phase_value_at_entrained_time = analytical_phase_difference(end);
%     light_profiles_index
    save_analytical_phase_difference_LD(:,light_profiles_index) = analytical_phase_difference;
    assignin('base','save_analytical_phase_difference_LD',save_analytical_phase_difference_LD)
    assignin('base','POE_ROE_Avg',[entrained_time,phase_value_at_entrained_time,average_POE])
    assignin('base','analytical_phase_difference',analytical_phase_difference)



function [ derivative ] = RK4_circadian_CIRC(y)
    theta_I = y(1,1);
    omega_I(i) = y(1,2);
    Normalized_theta_Z(i) = mod(zeitgeber_clock(i),2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);

    % Light profile
    if abs(Normalized_theta_Z(i) - light_middle) < epsilon 
        if mod(i,light_period/T_interval) > light_on/T_interval && mod(i,light_period/T_interval) < light_off / T_interval
            light(i) = C_ZG;
        else
            light(i) = 0;
        end
    else
        light(i) = 0;
    end
    
    % dCiRC function
    if Normalized_theta_I(i) >= CIRC_shift && Normalized_theta_I(i) <= pi + CIRC_shift
        CIRC(i) = sin(Normalized_theta_I(i) - CIRC_shift) + s * sin(2 * (Normalized_theta_I(i) - CIRC_shift));
        if CIRC(i) < 0
            CIRC(i) = 0;
        end
        if a < 1
            CIRC(i) = CIRC(i) * a;
        end
    else
        CIRC(i) = -sin(2 * pi - (Normalized_theta_I(i) - CIRC_shift)) -...
               s * sin(2 * pi - (2 * (Normalized_theta_I(i) - CIRC_shift)));
        if CIRC(i) > 0
            CIRC(i) = 0;
        end
        if a > 1
            CIRC(i) = CIRC(i) / a;
        end
    end

    if s < 1
        max_1 = -0.5669 * s^3 + 1.1431 * s^2 + 0.1703 * s + 0.9963;
        CIRC(i) = CIRC(i) / max_1;
    else
        max_2 = 0.0029 * s^2 + 0.9738 * s + 0.7783;
        CIRC(i) = CIRC(i) / max_2;
    end
    
    % dCiRC Model ODEs
        dtheta_I = omega_I(i);
        domega_I(i) = light(i) * CIRC(i) + elastic_factor * (omega_DD - omega_I(i));
        derivative = [dtheta_I domega_I(i)];
end
end

function CIRC_LD_Cost_for_plot_best_fitting_DD(optimized_x,organism_index,T_interval,number_of_timesteps)
global T_final light_profiles_DD light_profiles_index zeit_t ...
       light_period light_off light_on save_analytical_phase_difference_DD

% preallocate
ALL_clocks                  = zeros(number_of_timesteps,2);
zeitgeber_clock             = zeros(number_of_timesteps,1);
Normalized_theta_Z          = zeros(number_of_timesteps,1);
Normalized_theta_I          = zeros(number_of_timesteps,1);
CIRC                        = zeros(number_of_timesteps,1);
light                       = zeros(number_of_timesteps,1);
omega_I                     = zeros(number_of_timesteps,1);
domega_I                    = zeros(number_of_timesteps,1);

% Parameters
tau_DD          = optimized_x(organism_index,1);
C_ZG            = optimized_x(organism_index,2);
elastic_factor  = optimized_x(organism_index,3);
s               = optimized_x(organism_index,4);
a               = optimized_x(organism_index,5);
phi_0           = optimized_x(organism_index,6);
zeit_shift      = optimized_x(organism_index,7);
omega_DD        = 2 * pi ./ tau_DD;
initial_theta_I = phi_0;
initial_omega_I = omega_DD;
CIRC_shift      = 0;

% light profiles
epsilon              = light_profiles_DD(light_profiles_index,1); % radius of entrainment window
light_middle         = light_profiles_DD(light_profiles_index,2); % center of entrainment window

%% Simulations
theta_I_initial = initial_theta_I; 
omega_I_initial = initial_omega_I; 
initial_conditions = [theta_I_initial, omega_I_initial];
ALL_clocks(1,:) = initial_conditions;

for i = 1:number_of_timesteps - 1
    zeitgeber_clock(i) = (i-1-(18-zeit_shift)/T_interval)*2*pi/(zeit_t/T_interval);
    k1= RK4_circadian_CIRC(ALL_clocks(i,:));
    k2_clocks = ALL_clocks(i,:) + 0.5*T_interval*k1;

    k2= RK4_circadian_CIRC(k2_clocks);
    k3_clocks = ALL_clocks(i,:) + 0.5*T_interval*k2;

    k3= RK4_circadian_CIRC(k3_clocks);
    k4_clocks = ALL_clocks(i,:) + T_interval*k3;

    k4= RK4_circadian_CIRC(k4_clocks);
    ALL_clocks(i+1,:) =  ALL_clocks(i,:) + T_interval*(k1 + 2*k2 + 2*k3 + k4) / 6.0;             
end
    analytical_zeitgeber_signal = sin(zeitgeber_clock);
    theta_I_signal = sin(ALL_clocks(:,1));

% find the time values corresponding to the peaks of the signal
    [~, n_time]= findpeaks(theta_I_signal,'MinPeakHeight',0.99);
    [~, analytical_z_time]= findpeaks(analytical_zeitgeber_signal,'MinPeakHeight',0.99);

    analytical_days = size(analytical_z_time,1);
    analytical_phase_difference = zeros(analytical_days,1);
    
    for j = 1:analytical_days
        [~, analytical_phase_index] = min(abs(analytical_z_time(j) - n_time));
        analytical_phase_difference(j) = (n_time(analytical_phase_index) - analytical_z_time(j)) .* T_interval;
    end
    
    entrained_time = T_final;
    average_POE = mean(analytical_phase_difference);
    phase_value_at_entrained_time = analytical_phase_difference(end);
%     light_profiles_index
    save_analytical_phase_difference_DD(:,light_profiles_index) = analytical_phase_difference;
    assignin('base','save_analytical_phase_difference_DD',save_analytical_phase_difference_DD)
    assignin('base','POE_ROE_Avg',[entrained_time,phase_value_at_entrained_time,average_POE])
    assignin('base','analytical_phase_difference',analytical_phase_difference)


        
function [ derivative ] = RK4_circadian_CIRC(y)
    theta_I = y(1,1);
    omega_I(i) = y(1,2);
    Normalized_theta_Z(i) = mod(zeitgeber_clock(i),2*pi);
    Normalized_theta_I(i) = mod(theta_I,2*pi);

    % Light profile
    if abs(Normalized_theta_Z(i) - light_middle) < epsilon 
        if mod(i,light_period/T_interval) > light_on/T_interval && mod(i,light_period/T_interval) < light_off / T_interval
            light(i) = C_ZG;
        else
            light(i) = 0;
        end
    else
        light(i) = 0;
    end
    
    % dCiRC function
    if Normalized_theta_I(i) >= CIRC_shift && Normalized_theta_I(i) <= pi + CIRC_shift
        CIRC(i) = sin(Normalized_theta_I(i) - CIRC_shift) + s * sin(2 * (Normalized_theta_I(i) - CIRC_shift));
        if CIRC(i) < 0
            CIRC(i) = 0;
        end
        if a < 1
            CIRC(i) = CIRC(i) * a;
        end
    else
        CIRC(i) = -sin(2 * pi - (Normalized_theta_I(i) - CIRC_shift)) -...
               s * sin(2 * pi - (2 * (Normalized_theta_I(i) - CIRC_shift)));
        if CIRC(i) > 0
            CIRC(i) = 0;
        end
        if a > 1
            CIRC(i) = CIRC(i) / a;
        end
    end

    if s < 1
        max_1 = -0.5669 * s^3 + 1.1431 * s^2 + 0.1703 * s + 0.9963;
        CIRC(i) = CIRC(i) / max_1;
    else
        max_2 = 0.0029 * s^2 + 0.9738 * s + 0.7783;
        CIRC(i) = CIRC(i) / max_2;
    end
    % dCiRC Model ODEs
        dtheta_I = omega_I(i);
        domega_I(i) = light(i) * CIRC(i) + elastic_factor * (omega_DD - omega_I(i));
        derivative = [dtheta_I domega_I(i)];
end
end